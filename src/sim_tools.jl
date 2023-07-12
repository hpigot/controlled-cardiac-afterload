using DelimitedFiles, Polynomials, Plots, ControlSystems

"[m^3/s] per [L/min]"
kϕ = 1 / (60 * 1000)
"[Pa] per [mmHg]"
kp = 133.3
"Pa to mmHg"
convert_p(p) = p / kp
"m*3/s to L/min"
convert_ϕ(ϕ) = ϕ / kϕ

"Identified parallel 4-element windkessel human-data model parameters from pigot2021identification, 10.1016/j.ifacol.2021"
function θwk4()
    Rp = 13.6 # mmHg/(L/min) Peripheral (systemic) resistance
    C = 0.0743 # L/mmHg Total arterial compliance
    Rc = 0.952 # mmHg/(L/min) characeristic aoritic resistance
    L = 0.0952 # mmHg*min/L Total arterial impedance
    return Rp, C, Rc, L
end

"Identified parallel 2-element windkessel human-data model parameters from pigot2021identification, 10.1016/j.ifacol.2021"
function θwk2()
    Rp = 13.6 # mmHg/(L/min) Peripheral (systemic) resistance
    C = 0.0996 # L/mmHg Total arterial compliance
    return Rp, C
end

"returns pressure [mmHg], with cardiac power w [W] and plunger position z [m]"
function get_p(z, w; plast=0.0)
    y = y0 * z / zmax # admittance [(L/min)/mmHg]
    w_ = w / (kp * kϕ) # [mmHg*L/min]
    p = Polynomial([-w_ / y, -ϕaux / y, 1]) |> roots

    if plast <= 0 # first p solution where w>0
        plast = maximum(p[findfirst(x -> x > 0, w)])
    end

    p[argmin(abs.(p .- plast))]
end

"Aortic flow based on pressure p [mmHg] and plunger position z [m]"
function get_ϕ(z, p)
    y0 * z * p / zmax - ϕaux
end

"Cardiac power [mmHg*L/min] based on pressure p [mmHg] and plunger position z [m]"
function get_w(z, p)
    y0 * z * p^2 / zmax - ϕaux * p
end

"plunger position z [m] that yeilds pressure p [mmHg] from flow ϕ [L/min]."
function get_z(p, ϕ)
    zmax * (ϕ + ϕaux) / (y0 * p)
end

"Given current plunger position `z`, pressure `pk`, and previous pressures and flows `pks` and `ϕks`, calculate the desired pressure, flow, and reslting plunger position that satisfy 4-element Windkessel dynamics with transfer function denominator and numerator coefficients `a` and `b` (defined globally)."
function simulate_wk4(z, pk, pks, ϕks)
    α = a[1] / b[1]
    β = (a[2] * pks[end] + a[3] * pks[end-1] - b[2] * ϕks[end] - b[3] * ϕks[end-1]) / b[1]
    wk = y0 * z * pk^2 / zmax - ϕaux * pk
    p = Polynomial([-wk, β, α]) |> roots
    p_nearest = p[argmin(abs.(p .- pks[end]))]
    ϕ = α * p_nearest + β
    z = get_z(p_nearest, ϕ)
    p_nearest, ϕ, z
end

"Given current plunger position `z`, pressure `pk`, and previous pressures and flows `pks` and `ϕks`, calculate the desired pressure, flow, and reslting plunger position that satisfy 2-element Windkessel dynamics with model paramters `C2` and `Rp2` (defined globally)."
function simulate_wk2(z, pk, pks, ϕks)
    τ = -h / (C2 * Rp2)
    p = Rp2 * (1 - exp(τ)) * ϕks[end] + exp(τ) * pks[end]
    wk = get_w(z, pk)
    ϕ = wk / p
    z = get_z(p, ϕ)
    p, ϕ, z
end

mutable struct AfterloadSimulation
    p::Vector{Float64}
    ϕ::Vector{Float64}
    z::Vector{Float64}
    AfterloadSimulation(length::Int) = new(Vector{Float64}(undef, length), Vector{Float64}(undef, length), Vector{Float64}(undef, length))
end

function simulate_afterload(w::Vector{Float64}, dynamics::Function, plims; z0=z0, zmin=zmin, ϕaux=ϕaux)
    measured = AfterloadSimulation(length(w))
    result = AfterloadSimulation(length(w))
    result_lim = AfterloadSimulation(length(w))

    # initial conditions
    initlen = 2
    p0 = get_p.(fill(z0, initlen), w[1:initlen], plast=500) # first two time steps
    !isnothing(findfirst(x -> (x < 0), p0)) ? (@error "p0 negative") : nothing
    ϕ0 = get_ϕ.(fill(z0, initlen), p0) # first two time steps

    # initialize simulation vectors
    measured.p[1:initlen] = p0
    measured.ϕ[1:initlen] = ϕ0
    result.p[1:initlen] = p0
    result.ϕ[1:initlen] = ϕ0
    result.z[1:initlen] = fill(z0, initlen)
    result_lim.p[1:initlen] = p0
    result_lim.ϕ[1:initlen] = ϕ0
    result_lim.z[1:initlen] = fill(z0, initlen)

    for k = 3:length(w)
        # get current pressure measurement, with z position held from previous actuation
        zk = result_lim.z[k-1]
        pk = get_p(zk, w[k], plast=result_lim.p[k-1]) # may cause complex roots.
        ϕk = get_ϕ(zk, pk)

        if !isreal(pk)
            @warn "solution failed at k=$(k), pk = $(pk), ϕk = $(ϕk), zk = $(zk)"
            zk *= 0.5
            pk = get_p(zk, w[k], plast=result.p[k-1])
        end

        measured.p[k] = pk
        measured.ϕ[k] = ϕk

        # controlled system, if dynamics are followed without clamping
        pref, ϕref, zref = dynamics(zk, pk, result.p[k-2:k-1], result.ϕ[k-2:k-1])

        result.p[k] = pref
        result.ϕ[k] = ϕref
        result.z[k] = zref

        # clamp pressure
        if k > 0
            wk = get_w(zk, measured.p[k])
            pref = clamp(pref, plims...)
            ϕref = wk / pref
            zref = get_z(pref, ϕref)
        end

        # assert min zref
        if zref < zmin
            @show "zref $(zref) clamped to $(zmin)"
            zref = zmin
        end

        # assert that the flow does not exceed ϕaux 
        if ϕref < -ϕaux
            @warn "k=$(k)ϕref < -ϕaux\t$(ϕref) < $(-ϕaux)"
            ϕaux = max(convert_ϕ(1.1 * -ϕref), 0)
        end

        result_lim.z[k] = zref
        result_lim.p[k] = pref
        result_lim.ϕ[k] = ϕref
    end

    result, result_lim, measured
end

"return repetition of `xbeat` equal to the length of xbeat repeated `nbeat` times, with the beat number `shiftat` shifted forward by a fraction `shiftby` of the beat period."
function shifted_beat(xbeat, nbeats, shiftat, shiftby)
    ibeat = length(xbeat)
    x = repeat(xbeat, nbeats)
    ishift = Int(round((shiftat - shiftby) * ibeat))
    x = vcat(x[1:ishift], repeat(xbeat, nbeats - shiftat + 1))
    x[1:nbeats*ibeat]
end

"return repetition of `xbeat` equal to the length of xbeat repeated `nbeat` times, with the beat number `removeat` removed, and written over by the end beat value. If `foh`, then written over accordnig to first order hold between the start and end beat values."
function missing_beat(xbeat, nbeats, removeat; foh=false)
    ibeat = length(xbeat)
    x = repeat(xbeat, nbeats)
    irm = Int(round(ibeat * (removeat - 1)))
    irm:irm+ibeat-1
    if foh
        @show x[irm:irm+ibeat-1] = [xbeat[1] + i * (xbeat[end] - xbeat[1]) / ibeat for i in 0:ibeat-1]
    else
        x[irm:irm+ibeat-1] .= xbeat[end]
    end
    x
end

function plot_all(t, w, result, result_lim, measured, reference; label="", plotextra=false)
    pp = plot(t, result_lim.p, label="p⁺ limited",ylabel="Pressure [mmHg]")
    plot!(t, result.p, label="p⁺ " * label)
    pϕ = plot(t, result_lim.ϕ, label="ϕ⁺ limited", ylabel="Flow [L/min]")
    plot!(t, result.ϕ, label="ϕ⁺ " * label)
    pz = plot(t, result_lim.z, label="z⁺ limited", ylabel="Position [m]")
    plot!(t, result.z, label="z⁺ " * label)
    pw = plot(t, w, label="w", ylabel="Power [W]", xlabel="Time [s]")

    if plotextra # include reference dynamics, simulated measured values and sampled input data plots
        plot!(pp, t, measured.p, label="p measured")
        plot!(pp, t, reference.p, label="p "*label*" ref")
        plot!(pp, t, p1, label="p human")
        plot!(pϕ, t, measured.ϕ, label="ϕ measured")
        plot!(pϕ, t, ϕ1, label="ϕ human")
        hline!(pϕ, [ϕaux], label="ϕaux")
    end

    plot(pp, pϕ, pz, pw, layout=(4, 1), size=(800, 800))
end

function plot_steady(istart, iend, t, w, result::AfterloadSimulation, result_lim::AfterloadSimulation, measured::AfterloadSimulation, reference::AfterloadSimulation; label="", plotextra=false)
    ts = t[istart:iend]
    pp = plot(ts, result_lim.p[istart:iend], label="p⁺ limited",ylabel="Pressure [mmHg]")
    plot!(ts, result.p[istart:iend], label="p⁺ " * label) 
    pϕ = plot(ts, result_lim.ϕ[istart:iend], label="ϕ⁺ limited", ylabel="Flow [L/min]")
    plot!(ts, result.ϕ[istart:iend], label="ϕ⁺ " * label)
    pz = plot(ts, result_lim.z[istart:iend], label="z⁺ limited", ylabel="Position [m]")
    plot!(ts, result.z[istart:iend], label="z⁺ " * label)
    pw = plot(ts, w[istart:iend], label="w", ylabel="Power [W]", xlabel="Time [s]")
    
    if plotextra # include reference dynamics, simulated measured values and sampled input data plots
        plot!(pp, ts, measured.p[istart:iend], label="p measured")
        plot!(pp, ts, reference.p[istart:iend], label="p "*label*" ref")
        plot!(pp, ts, p1[istart:iend], label="p human")
        plot!(pϕ, ts, measured.ϕ[istart:iend], label="ϕ measured")
        plot!(pϕ, ts, ϕ1[istart:iend], label="ϕ human")
        hline!(pϕ, [ϕaux], label="ϕaux")
    end

    plot(pp, pϕ, pz, pw, layout=(4, 1), size=(800, 800))
end

function plot_residuals(istart, iend, t, results, reference)
    ϵp = reference.p[istart:iend].-results.p[istart:iend]
    ϵϕ = reference.ϕ[istart:iend].-results.ϕ[istart:iend]
    ϵw = ((reference.ϕ[istart:iend].*reference.p[istart:iend]).-(results.ϕ[istart:iend].*results.p[istart:iend])) * kp * kϕ

    pp = plot(t[istart:iend],ϵp,title="max steady-state p error $(printn(maximum(abs.(ϵp))))", ylabel="Pressure [mmHg]",label="ϵp");
    pϕ = plot(t[istart:iend],ϵϕ,title="max steady-state ϕ error $(printn(maximum(abs.(ϵϕ))))", ylabel="Flow [L/min]",label="ϵϕ");
    pw = plot(t[istart:iend],ϵw,title="max steady-state w error $(printn(maximum(abs.(ϵw))))", ylabel="Power [W]",label="ϵw");

    plot(pp, pϕ, pw, layout=(3, 1), size=(600, 800))
end

"write matrix d to file `filename.csv`, with `header` and `comment`."
function csvwrite(filename::String, d::Matrix; comment::String="", header=["x$(n)" for n in 1:size(d)[2]])
    h = string([s * "," for s in header]...)[1:end-1]
    s = isempty(comment) ? h * "\n" : "# " * comment * "\n" * h * "\n"
    open(filename * ".csv", "w") do f
        write(f, s)
        writedlm(f, d, ',')
    end
end

"round number to 4 sigdigits for printing"
printn(n) = round(n, sigdigits=4)