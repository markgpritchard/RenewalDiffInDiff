
using DrWatson 
@quickactivate :RenewalDiffInDiff

using Distributions
using Random
using RenewalDiD 

## Simulation 1:
# 3 large populations 
# 0.1% initially exposed
# constant transmission, equal for all groups, R0
# constant detection, 80%
# no effective intervention; placebo intervention in one group at time 50
sim1 = let 
    rng = Xoshiro(1)
    Ns = [round(Int, rand(rng, Uniform(30_000_000, 40_000_000))) for _ in 1:2]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:2]
    eta = 0.2
    sigma = 0.5
    phi = 0.8
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta2_a = x -> beta2ratio * beta1(x)
    beta2_b = x -> beta2_a(x) * (x < 50 ? 1 : 0.8)
    betas_a = [beta1, beta2_a]
    betas_b = [beta1, beta2_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=50,
    )
    sima = packsimulations(rng, 100, s1a, s2a; id="sim1a", minvalue=0.5, sampletime=10,)
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=50,
    )
    simb = packsimulations(rng, 100, s1b, s2b; id="sim1b", minvalue=0.5, sampletime=10,)
    @ntuple sima simb
end
safesave(simulationdir("sim1.jld2"), Dict("sim" => sim1))

## Simulation 2:
# as in simulation 1 but intervention in one group at time 50 reduces transmission by 20%
sim2 = let 
    rng = Xoshiro(2)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = repeat([2 * mu]; inner=3)
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim2", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim2.jld2"), Dict("sim" => sim2))

## Simulation 3:
# as in simulation 1 but 10% variation in transmission between groups
sim3 = let 
    rng = Xoshiro(3)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas = [2 * mu * rand(rng, Uniform(0.9, 1.1)) for _ in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim3", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim3.jld2"), Dict("sim" => sim3))

## Simulation 4:
# as in simulation 3 but intervention in one group at time 50 reduces transmission by 20%
sim4 = let 
    rng = Xoshiro(4)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = [2 * mu * rand(rng, Uniform(0.9, 1.1)) for _ in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim4", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim4.jld2"), Dict("sim" => sim4))

## Simulation 5:
# as in simulation 1 but 25% variation in transmission between groups
sim5 = let 
    rng = Xoshiro(5)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim5", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim5.jld2"), Dict("sim" => sim5))

## Simulation 6:
# as in simulation 5 but intervention in one group at time 50 reduces transmission by 20%
sim6 = let 
    rng = Xoshiro(6)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim6", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim6.jld2"), Dict("sim" => sim6))

## Simulation 7:
# as in simulation 5 but 20% seasonal forcing of transmission
sim7 = let 
    rng = Xoshiro(7)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betas = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim7", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim7.jld2"), Dict("sim" => sim7))

## Simulation 8:
# as in simulation 7 but intervention in one group at time 50 reduces transmission by 20%
sim8 = let 
    rng = Xoshiro(8)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), t -> betasc[2](t), t -> betasc[3](t)]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim8", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim8.jld2"), Dict("sim" => sim8))


## Simulation 9:
# as in simulation 8 but intervention in two groups at time 50 reduces transmission by 20%
sim9 = let 
    rng = Xoshiro(9)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 50 ? 1.0 : 0.8) * betasc[2](t), 
        t -> betasc[3](t)
    ]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=50,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim9", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim9.jld2"), Dict("sim" => sim9))

## Simulation 10:
# as in simulation 9 but intervention in two groups at different times reduces transmission by 20%
sim10 = let 
    rng = Xoshiro(10)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * betasc[2](t), 
        t -> betasc[3](t)
    ]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=25,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim10", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim10.jld2"), Dict("sim" => sim10))

## Simulation 11:
# as in simulation 10 but intervention of interest has no effect and an alternative 
# intervention reduces transmission by 20%
sim11 = let 
    rng = Xoshiro(11)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim11", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim11.jld2"), Dict("sim" => sim11))

## Simulation 12:
# as in simulation 11 but intervention of interest reduces transmission by 20%
sim12 = let 
    rng = Xoshiro(12)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim12", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim12.jld2"), Dict("sim" => sim12))

## Simulation 13:
# as in simulation 11 but alternative intervention increases transmission by 20%
sim13 = let 
    rng = Xoshiro(13)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 1.2) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 1.2) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim13", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim13.jld2"), Dict("sim" => sim13))

## Simulation 14:
# as in simulation 13 but intervention of interest reduces transmission by 20%
sim14 = let 
    rng = Xoshiro(14)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 1.2) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 1.2) * betasc[3](t)
    ]

    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim14", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim14.jld2"), Dict("sim" => sim14))

## Simulation 15:
# as in simulation 11 but lower proportion detected
sim15 = let 
    rng = Xoshiro(15)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.4
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim15", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim15.jld2"), Dict("sim" => sim15))

## Simulation 16:
# as in simulation 15 but intervention of interest reduces transmission by 20%
sim16 = let 
    rng = Xoshiro(16)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.4
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim16", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim16.jld2"), Dict("sim" => sim16))

## Simulation 17:
# as in simulation 15 but proportion detected changes over time, consistently for all groups
sim17 = let 
    rng = Xoshiro(17)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = t -> 0.4 + 0.004 * t
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim17", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim17.jld2"), Dict("sim" => sim17))

## Simulation 18:
# as in simulation 17 but intervention of interest reduces transmission by 20%
sim18 = let 
    rng = Xoshiro(18)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = t -> 0.4 + 0.004 * t
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim18", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim18.jld2"), Dict("sim" => sim18))

## Simulation 19:
# as in simulation 15 but alternative intervention increases proportion detected
sim19 = let 
    rng = Xoshiro(19)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> 0.4, 
        t -> (t < 60 ? 0.4 : 0.8), 
        t -> (t < 35 ? 0.4 : 0.8)
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim19", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim19.jld2"), Dict("sim" => sim19))

## Simulation 20:
# as in simulation 19 but intervention of interest reduces transmission by 20%
sim20 = let 
    rng = Xoshiro(20)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> 0.4, 
        t -> (t < 60 ? 0.4 : 0.8), 
        t -> (t < 35 ? 0.4 : 0.8)
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim20", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim20.jld2"), Dict("sim" => sim20))

## Simulation 21:
# as in simulation 15 but intervention of interest increases proportion detected
sim21 = let 
    rng = Xoshiro(21)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> (t < 50 ? 0.4 : 0.8), 
        t -> (t < 24 ? 0.4 : 0.8), 
        t -> 0.4
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim21", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim21.jld2"), Dict("sim" => sim21))

## Simulation 22:
# as in simulation 21 but intervention of interest reduces transmission by 20%
sim22 = let 
    rng = Xoshiro(22)
    u0s = [simu0(rng, largepop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> (t < 50 ? 0.4 : 0.8), 
        t -> (t < 24 ? 0.4 : 0.8), 
        t -> 0.4
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim22", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim22.jld2"), Dict("sim" => sim22))

## Simulation 23:
# as in simulation 1 but with small populations
sim23 = let 
    rng = Xoshiro(23)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas = repeat([2 * mu]; inner=3)
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim23", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim23.jld2"), Dict("sim" => sim23))

## Simulation 24:
# as in simulation 2 but with small populations
sim24 = let 
    rng = Xoshiro(24)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = repeat([2 * mu]; inner=3)
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim24", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim24.jld2"), Dict("sim" => sim24))

## Simulation 25:
# as in simulation 3 but with small populations
sim25 = let 
    rng = Xoshiro(25)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas = [2 * mu * rand(rng, Uniform(0.9, 1.1)) for _ in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim25", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim25.jld2"), Dict("sim" => sim25))

## Simulation 26:
# as in simulation 4 but with small populations
sim26 = let 
    rng = Xoshiro(26)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = [2 * mu * rand(rng, Uniform(0.9, 1.1)) for _ in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim26", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim26.jld2"), Dict("sim" => sim26))

## Simulation 27:
# as in simulation 5 but with small populations
sim27 = let 
    rng = Xoshiro(27)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim27", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim27.jld2"), Dict("sim27" => sim5))

## Simulation 28:
# as in simulation 6 but with small populations
sim28 = let 
    rng = Xoshiro(28)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betasc = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1], t -> betasc[2], t -> betasc[3]]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim28", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim28.jld2"), Dict("sim" => sim28))

## Simulation 29:
# as in simulation 7 but with small populations
sim29 = let 
    rng = Xoshiro(29)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betas = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim29", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim29.jld2"), Dict("sim" => sim29))

## Simulation 30:
# as in simulation 8 but with small populations
sim30 = let 
    rng = Xoshiro(30)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), t -> betasc[2](t), t -> betasc[3](t)]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=nothing,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim30", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim30.jld2"), Dict("sim" => sim30))

## Simulation 31:
# as in simulation 9 but with small populations
sim31 = let 
    rng = Xoshiro(31)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 50 ? 1.0 : 0.8) * betasc[2](t), 
        t -> betasc[3](t)
    ]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=50,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim31", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim31.jld2"), Dict("sim" => sim31))

## Simulation 32:
# as in simulation 10 but with small populations
sim32 = let 
    rng = Xoshiro(32)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * betasc[2](t), 
        t -> betasc[3](t)
    ]
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=25,
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim32", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim32.jld2"), Dict("sim" => sim32))

## Simulation 33:
# as in simulation 11 but with small populations
sim33 = let 
    rng = Xoshiro(33)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim33", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim33.jld2"), Dict("sim" => sim33))

## Simulation 34:
# as in simulation 12 but with small populations
sim34 = let 
    rng = Xoshiro(34)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim34", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim34.jld2"), Dict("sim" => sim34))

## Simulation 35:
# as in simulation 13 but with small populations
sim35 = let 
    rng = Xoshiro(35)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 1.2) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 1.2) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim35", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim35.jld2"), Dict("sim" => sim35))

## Simulation 36:
# as in simulation 14 but with small populations
sim36 = let 
    rng = Xoshiro(36)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.8
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 1.2) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 1.2) * betasc[3](t)
    ]

    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim36", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim36.jld2"), Dict("sim" => sim36))

## Simulation 37:
# as in simulation 15 but with small populations
sim37 = let 
    rng = Xoshiro(37)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.4
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim37", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim37.jld2"), Dict("sim" => sim37))

## Simulation 38:
# as in simulation 16 but with small populations
sim38 = let 
    rng = Xoshiro(38)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = 0.4
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim38", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim38.jld2"), Dict("sim" => sim38))

## Simulation 39:
# as in simulation 17 but with small populations
sim39 = let 
    rng = Xoshiro(39)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = t -> 0.4 + 0.004 * t
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim39", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim39.jld2"), Dict("sim" => sim39))

## Simulation 40:
# as in simulation 18 but with small populations
sim40 = let 
    rng = Xoshiro(40)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psi = t -> 0.4 + 0.004 * t
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi, kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi, kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi, kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim40", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim40.jld2"), Dict("sim" => sim40))

## Simulation 41:
# as in simulation 19 but with small populations
sim41 = let 
    rng = Xoshiro(41)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> 0.4, 
        t -> (t < 60 ? 0.4 : 0.8), 
        t -> (t < 35 ? 0.4 : 0.8)
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim41", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim41.jld2"), Dict("sim" => sim41))

## Simulation 42:
# as in simulation 20 but with small populations
sim42 = let 
    rng = Xoshiro(42)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> 0.4, 
        t -> (t < 60 ? 0.4 : 0.8), 
        t -> (t < 35 ? 0.4 : 0.8)
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> (t < 50 ? 1.0 : 0.8) * betasc[1](t), 
        t -> (t < 25 ? 1.0 : 0.8) * (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim42", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim42.jld2"), Dict("sim" => sim42))

## Simulation 43:
# as in simulation 21 but with small populations
sim43 = let 
    rng = Xoshiro(43)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> (t < 50 ? 0.4 : 0.8), 
        t -> (t < 24 ? 0.4 : 0.8), 
        t -> 0.4
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim43", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim43.jld2"), Dict("sim" => sim43))

## Simulation 44:
# as in simulation 22 but with small populations
sim44 = let 
    rng = Xoshiro(44)
    u0s = [simu0(rng, smallpop, 0.001) for _ in 1:3]
    mu = 0.2
    kappa = 0.5
    delta = 0.3
    psis = [
        t -> (t < 50 ? 0.4 : 0.8), 
        t -> (t < 24 ? 0.4 : 0.8), 
        t -> 0.4
    ]
    betas_const = [2 * mu * rand(rng, Uniform(0.75, 1.25)) for _ in 1:3]
    betasc = [t -> betas_const[i] * (1 + 0.2 * cos(2pi * (t - 20) / 365)) for i in 1:3]
    betas = [
        t -> betasc[1](t), 
        t -> (t < 60 ? 1.0 : 0.8) * betasc[2](t), 
        t -> (t < 35 ? 1.0 : 0.8) * betasc[3](t)
    ]
    
    s1 = packsimulationtuple( ; 
        u0=u0s[1], beta=betas[1], mu, delta, psi=psis[1], kappa, intervention=[50, nothing],
    )
    s2 = packsimulationtuple( ; 
        u0=u0s[2], beta=betas[2], mu, delta, psi=psis[2], kappa, intervention=[25, 60],
    )
    s3 = packsimulationtuple( ; 
        u0=u0s[3], beta=betas[3], mu, delta, psi=psis[3], kappa, intervention=[nothing, 35],
    )
    packsimulations(rng, 100, s1, s2, s3; id="sim44", minvalue=0.1, sampletime=10,)
end
safesave(simulationdir("sim44.jld2"), Dict("sim" => sim44))
