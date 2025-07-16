
using CairoMakie 
using Random


## simple stochastic small steps 

ts = 0:1:1000

ys1 = zeros(length(ts), 20)
ys2 = ones(length(ts), 20)

for i ∈ 1:20, t ∈ ts 
    if t == 0 
        ys1[t+1, i] = 1.0 
    elseif rand() < 0.1 
        ys1[t+1, i] = max(0.0, ys1[t, i] - 0.01)
    else
        ys1[t+1, i] = ys1[t, i]
    end 
end

for i ∈ 1:20, t ∈ ts 
    t == 0 && continue 
    if rand() < 0.01 
        res = 0.01
        for ti ∈ t:1000 
            ys2[ti+1, i] = max(0.0, ys2[ti+1, i] - res)
            res = 0.01 * exp(log(1 / 0.01) * (1 - exp(-0.002 * (ti - t))))
        end
    end
end

fig = Figure()

ax1 = Axis(fig[1, 1])
for i ∈ 1:19 
    lines!(ax1, ts, ys1[:, i]; color=:lightgray, linewidth=1,)
end
lines!(ax1, ts, ys1[:, 20]; color=:black, linewidth=1,)

ax2 = Axis(fig[1, 2])
for i ∈ 1:19 
    lines!(ax2, ts, ys2[:, i]; color=:lightgray, linewidth=1,)
end
lines!(ax2, ts, ys2[:, 20]; color=:black, linewidth=1,)

Label(fig[1, 0], "Proportion susceptible"; fontsize=11.84, rotation=-π/2, tellheight=false)
Label(fig[2, 1:2], "time"; fontsize=11.84, tellwidth=false)

fig

save("fi1.pdf", fig)
