
using DrWatson
@quickactivate :RenewalDiffInDiff

using CairoMakie
using Distributions
using Random

xs = 1:1:100 

rng = Xoshiro(1729)

poissonys = [[rand(rng, Poisson(x)) for x in xs] for _ in 1:10]
normalys = [[(z = rand(rng, Normal(x, sqrt(x))); max(z, 0)) for x in xs] for _ in 1:10]
normalys2 = [[(z = rand(rng, Normal(0, 1)); x + z * sqrt(x)) for x in xs] for _ in 1:10]

fig = Figure() 
axs = [Axis(fig[1, i]) for i in 1:3]
for k in 1:10
    lines!(axs[1], xs, poissonys[k])
    lines!(axs[2], xs, normalys[k])
    lines!(axs[3], xs, normalys2[k])
end
linkaxes!(axs...)

fig

##

phis = [0.1, 0.5, 0.9]

binomys = [[rand(rng, Binomial(x, p)) for x in xs, p in phis] for _ in 1:10]
normalys = [[(z = rand(rng, Normal(x * p, x * p * (1 - p))); max(z, 0)) for x in xs, p in phis] for _ in 1:10]
normaly2s = [[(z = rand(rng, Normal(x * p, sqrt(x * p * (1 - p)))); max(z, 0)) for x in xs, p in phis] for _ in 1:10]

fig = Figure() 
axs = [Axis(fig[j, i]) for j in 1:3, i in 1:3]
for k in 1:10, j in 1:3
    lines!(axs[j, 1], xs, binomys[k][:, j])
    lines!(axs[j, 2], xs, normalys[k][:, j])
    lines!(axs[j, 3], xs, normaly2s[k][:, j])
end
linkaxes!(axs...)

fig

