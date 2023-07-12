include("sim_tools.jl")

plotextra = false # plot the results as displayed in manuscript
savefigures = false
savecsv = true
fileprefix = "csv_out/"

# Windkessel fit parameters from pigot2021identification, 10.1016/j.ifacol.2021.10.307
(Rp, C, Rc, L) = θwk4()
(Rp2, C2) = θwk2() # used in simulate_wk2()

# Human data, digitized from stergiopulos1999total,
# 10.1152/ajpheart.1999.276.1.H81 Figure 4, Type A, A and B
# sampled at 200Hz with Webplotdigitizer: Version 4.5
prefix = "data/"
p_sampled = readdlm(prefix * "pressure.csv", ',', Float64, '\n', comments=true, comment_char='#')[:, 2] # [mmHg]
ϕ_sampled = convert_ϕ.((1 / 1e6) .* readdlm(prefix * "flow.csv", ',', Float64, '\n', comments=true, comment_char='#')[:, 2]) # [L/min]

h = 0.005
t_sampled = collect(h .* (0:length(p_sampled)-1))
Tbeat = t_sampled[end]
ibeat = length(t_sampled)
totalbeats = 20
t1 = 0:h:h*totalbeats*ibeat-h

# repeated measured human data
ϕ1 = repeat(ϕ_sampled, totalbeats) # [L/min]
p1 = repeat(p_sampled, totalbeats) # [mmHg]
# @assert convert_ϕ(minimum(ϕc.(t1))) > convert_ϕ(-ϕaux) "$(convert_ϕ(minimum(ϕc.(t1)))) !> $(convert_ϕ(-ϕaux)), increase |ϕaux|"
w1 = ϕ1 .* p1 * kϕ * kp # [W]
# plot(t1, [ϕ1 p1 w1], labels=["ϕ1" "p1" "w1"], title="Continuous-time simulated WK4 data")

# discrete WK4 model
s = tf("s")
h = 0.005
Gc = Rc + Rp / (1 + s * C * Rp) - Rc / (1 + s * L / Rc) # Windkessel model
Gd = c2d(Gc, h) # [mmHg/L/min]
b = num(Gd)[1] # used in simulate_wk4()
a = den(Gd)[1] # used in simulate_wk4()

# pressure based on discrete WK4 model and measured flow as input
pwk4p, _, x, _ = lsim(Gd, (ϕ1)', collect(t1), x0=[0, 0], method=:zoh) # with transient
pwk4p, _, x, _ = lsim(Gd, (ϕ1)', collect(t1), x0=x[:, end], method=:zoh) # steady state
pwk4 = pwk4p'[:] # pressure from discrete WK4 [mmHg]
plot(pwk4, title="4 element Windkessel model", ylabel="Pressure [mmHg]", label="none")


# controlled cardiac afterload model parameters
γ = 12 * (Rp)^-1 / 0.03 # inversly proportional to plunger displacement z, set heuristically to give z < 0.03 with ϕaux at 21 L/min.
zmax = 0.03 # [m]
y0 = γ * zmax # [L/min/mmHg]
z0 = 0.5 * zmax # [m]
zmin = 1e-5 # [m]
ϕaux = 21 # [m^3/s] auxillary flow, note that -6.5 L/min is the minimum from the input signal; the larger value here avoids invalid solutions caused by w at the next time step using the previously calculated z

# simulate afterload controlled to WK4 dynamics
w4 = pwk4 .* ϕ1 * kp * kϕ # power from human flow measurments and resulting discrete wk4 pressures [W]
plims4 = (80, 105);
r4, rl4, m4 = simulate_afterload(w4, simulate_wk4, plims4)

ref4 = AfterloadSimulation(length(w4)) # The target values
ref4.p = @view pwk4[:]
ref4.ϕ = @view ϕ1[:]

plot_all(t1, w4, r4, rl4, m4, ref4, label="wk4", plotextra=plotextra) # including initial transient
nbeat = 2
istart = length(t1) - ibeat * nbeat - 1
iend = length(t1)
plot_steady(istart, iend, t1, w4, r4, rl4, m4, ref4, label="wk4", plotextra=plotextra)
savefigures ? savefig("sim_wk4.png") : nothing

plot_residuals(istart, iend, t1, r4, ref4)

# save wk4 sim data
if savecsv
    f4 = "sim_wk4"
    c4 = "4-element windkessel model, (Rp, C, Rc, L) = $(θwk4()).
    # auxiliary flow = $(convert_ϕ(ϕaux)|>printn) [L/min], zmax = $(printn(zmax)), y0 = $(printn(y0)) [L/min/mmHg], pressure limits $(plims4) [mmHg], h = $(h) [s]
    # t, time [s]; p, simulated aortic pressure [mmHg]; plim, simulated aortic pressure with limits [mmHg], f: simulated aortic flow [L/min]; flim: simulated aortic flow with limits [L/min]; z: afterload plunger position [m], zlim: afterload plunger position with limits [m], w: cardiac power [W]"
    h4 = ["t", "p", "plim", "f", "flim", "z", "zlim", "w"]
    t4 = round.(t1[istart:iend] .- t1[istart], digits=6)
    t4[1] = 0 # avoid -0.0
    d4 = hcat(t4, r4.p[istart:iend], rl4.p[istart:iend], r4.ϕ[istart:iend], rl4.ϕ[istart:iend], r4.z[istart:iend], rl4.z[istart:iend], w4[istart:iend])
    csvwrite(fileprefix * f4, d4, comment=c4, header=h4)
end

# discrete WK2 model
Gc2 = Rp2 / (1 + s * C2 * Rp2)
Gd2 = c2d(Gc2, h) # [mmHg/L/min]

pwk2p, _, x, _ = lsim(Gd2, (ϕ1)', collect(t1), x0=[0], method=:zoh) # [mmHg]
pwk2p, _, _, _ = lsim(Gd2, (ϕ1)', collect(t1), x0=x[:, end], method=:zoh) # [mmHg]
pwk2 = pwk2p'[:] # pressure from discrete WK4

# simulate afterload controlled to WK2 dynamics
w2 = pwk2 .* ϕ1 * kp * kϕ
plims2 = (80, 105);
r2, rl2, m2 = simulate_afterload(w2, simulate_wk2, plims2)

ref2 = AfterloadSimulation(length(w2))
ref2.p = @view pwk2[:]
ref2.ϕ = @view ϕ1[:]

plot_all(t1, w2, r2, rl2, m2, ref2, label="wk2", plotextra=plotextra) # including initial transient
istart = length(t1) - ibeat * nbeat - 1
iend = length(t1)
plot_steady(istart, iend, t1, w2, r2, rl2, m2, ref2, label="wk2", plotextra=plotextra)

savefigures ? savefig("sim_wk2.png") : nothing

plot_residuals(istart, iend, t1, r2, ref2)

if savecsv
    f2 = "sim_wk2"
    c2 = "2-element windkessel model, (Rp, C) = $(θwk2()).
    # auxiliary flow = $(convert_ϕ(ϕaux)|>printn) [L/min], zmax = $(printn(zmax)), y0 = $(printn(y0)) [L/min/mmHg], pressure limits $(plims2) [mmHg], h = $(h) [s]
    # t, time [s]; p, simulated aortic pressure [mmHg]; plim, simulated aortic pressure with limits [mmHg], f: simulated aortic flow [L/min]; flim: simulated aortic flow with limits [L/min]; z: afterload plunger position [m], zlim: afterload plunger position with limits [m], w: cardiac power [W]"
    h2 = ["t", "p", "plim", "f", "flim", "z", "zlim", "w"]
    t2 = round.(t1[istart:iend] .- t1[istart], digits=6)
    t2[1] = 0 # avoid -0.0
    d2 = hcat(t2, r2.p[istart:iend], rl2.p[istart:iend], r2.ϕ[istart:iend], rl2.ϕ[istart:iend], r2.z[istart:iend], rl2.z[istart:iend], w4[istart:iend])
    csvwrite(fileprefix * f2, d2, comment=c2, header=h2)
end

# initialization behavior
tinit = t1[1:3*ibeat]
p4init = plot(tinit, pwk4[1:3*ibeat], label="pwk4", ylabel="p [mmHg]", xlabel="t [s]");
plot!(tinit, r4.p[1:3*ibeat], label="p");
p2init = plot(tinit, pwk2[1:3*ibeat], label="pwk2", ylabel="p [mmHg]", title="Initial Transient Behavior");
plot!(tinit, r2.p[1:3*ibeat], label="p");
plot(p2init, p4init, layout=(2, 1))

if savecsv
    finit = "sim_init"
    cinit = "2-element windkessel model, (Rp, C) = $(θwk2()).
    # 4-element windkessel model, (Rp, C, Rc, L) = $(θwk4()).
    # auxiliary flow = $(convert_ϕ(ϕaux)|>printn) [L/min], zmax = $(printn(zmax)), y0 = $(printn(y0)) [L/min/mmHg], h = $(h) [s]
    # t, time [s]; pwk2, 2-element Windkessel aortic pressure [mmHg]; pswk2, simulated afterload aortic pressure with wk2 dynamics [mmHg]; pwk4, 4-element Windkessel aortic pressure [mmHg]; pswk4, simulated afterload aortic pressure with wk4 dynamics [mmHg]."
    hinit = ["t", "pwk2", "pswk2", "pwk4", "pswk4"]
    dinit = hcat(tinit, pwk2[1:3*ibeat], r2.p[1:3*ibeat], pwk4[1:3*ibeat], r4.p[1:3*ibeat])
    csvwrite(fileprefix * finit, dinit, comment=cinit, header=hinit)
end

# arrythmia 1: double beat
wd = shifted_beat(ϕ_sampled .* p_sampled * kϕ * kp, 20, 18, 0.5) # [W]
plimsd = (50, 120);
rd, rld, md = simulate_afterload(wd, simulate_wk4, plimsd)

plot_all(t1, wd, rd, rld, md, ref4, label="double", plotextra=plotextra) # including initial transient
nbeat = 4
istart = length(t1) - ibeat * nbeat - 1
iend = length(t1)
plot_steady(istart, iend, t1, wd, rd, rld, md, ref4, label="double", plotextra=plotextra)
savefigures ? savefig("sim_double.png") : nothing

# power residuals
ϵwd = wd[istart:iend] .- (rd.ϕ[istart:iend] .* rd.p[istart:iend]) * kp * kϕ
plot(t1[istart:iend], ϵwd, title="max steady-state w error $(printn(maximum(abs.(ϵwd))))", ylabel="Power [W]", label="ϵw double")

# save double beat sim data
if savecsv
    fd = "sim_double"
    cd = "4-element windkessel model, (Rp, C, Rc, L) = $(θwk4()).
    # auxiliary flow = $(convert_ϕ(ϕaux)|>printn) [L/min], zmax = $(printn(zmax)), y0 = $(printn(y0)) [L/min/mmHg], pressure limits $(plims4) [mmHg], h = $(h) [s]
    # t, time [s]; p, simulated aortic pressure [mmHg]; plim, simulated aortic pressure with limits [mmHg], f: simulated aortic flow [L/min]; flim: simulated aortic flow with limits [L/min]; z: afterload plunger position [m], zlim: afterload plunger position with limits [m], w: cardiac power [W]"
    hd = ["t", "p", "plim", "f", "flim", "z", "zlim", "w"]
    td = round.(t1[istart:iend] .- t1[istart], digits=6)
    td[1] = 0 # avoid -0.0
    dd = hcat(td, rd.p[istart:iend], rld.p[istart:iend], rd.ϕ[istart:iend], rld.ϕ[istart:iend], rd.z[istart:iend], rld.z[istart:iend], wd[istart:iend])
    csvwrite(fileprefix * fd, dd, comment=cd, header=hd)
end

# arrythmia 2: missing beat
wm = missing_beat(ϕ_sampled .* p_sampled * kϕ * kp, 20, 18)
plimsm = (50, 120);
rm, rlm, mm = simulate_afterload(wm, simulate_wk4, plimsm)

plot_all(t1, wm, rm, rlm, mm, ref4, label="missing", plotextra=plotextra) # including initial transient
plot_steady(istart, iend, t1, wm, rm, rlm, mm, ref4, label="missing", plotextra=plotextra)
savefigures ? savefig("sim_missing.png") : nothing

# power residuals
ϵwm = wm[istart:iend] .- (rm.ϕ[istart:iend] .* rm.p[istart:iend]) * kp * kϕ
plot(t1[istart:iend], ϵwm, title="max steady-state w error $(printn(maximum(abs.(ϵwm))))", ylabel="Power [W]", label="ϵw missing")

# save missing beat sim data
if savecsv
    fm = "sim_missing"
    cm = "4-element windkessel model, (Rp, C, Rc, L) = $(θwk4()).
    # auxiliary flow = $(convert_ϕ(ϕaux)|>printn) [L/min], zmax = $(printn(zmax)), y0 = $(printn(y0)) [L/min/mmHg], pressure limits $(plims4) [mmHg], h = $(h) [s]
    # t, time [s]; p, simulated aortic pressure [mmHg]; plim, simulated aortic pressure with limits [mmHg], f: simulated aortic flow [L/min]; flim: simulated aortic flow with limits [L/min]; z: afterload plunger position [m], zlim: afterload plunger position with limits [m], w: cardiac power [W]"
    hm = ["t", "p", "plim", "f", "flim", "z", "zlim", "w"]
    tm = round.(t1[istart:iend] .- t1[istart], digits=6)
    tm[1] = 0 # avoid -0.0
    dm = hcat(td, rm.p[istart:iend], rlm.p[istart:iend], rm.ϕ[istart:iend], rlm.ϕ[istart:iend], rm.z[istart:iend], rld.z[istart:iend], wm[istart:iend])
    csvwrite(fileprefix * fm, dm, comment=cm, header=hm)
end

# worst-case lower bound on axiliary flow value
ϕauxmax4 = 2 * sqrt(y0 * -minimum(w4 / (kp * kϕ)))
ϕauxmax2 = 2 * sqrt(y0 * -minimum(w2 / (kp * kϕ)))