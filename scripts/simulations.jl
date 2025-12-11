
using DrWatson 
@quickactivate :RenewalDiffInDiff

using Distributions
using Random
using RenewalDiD 

## Simulation 1:
# 3 large populations 
# 0.01% initially exposed
# constant detection, 80%
sim1 = let 
    rng = Xoshiro(1)
    rnga = Xoshiro(1)
    rngb = Xoshiro(1)
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
    sima = packsimulations(rnga, 100, s1a, s2a; id="sim1a", minvalue=0.5, sampletime=10,)
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=50,
    )
    simb = packsimulations(rngb, 100, s1b, s2b; id="sim1b", minvalue=0.5, sampletime=10,)
    @ntuple sima simb
end
safesave(simulationdir("sim1.jld2"), Dict("sim" => sim1))

## Simulation 2:
# as in simulation 1 but smaller populations
sim2 = let 
    rng = Xoshiro(2)
    rnga = Xoshiro(2)
    rngb = Xoshiro(2)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:2]
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
    sima = packsimulations(rnga, 100, s1a, s2a; id="sim2a", minvalue=0.5, sampletime=10,)
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=50,
    )
    simb = packsimulations(rngb, 100, s1b, s2b; id="sim2b", minvalue=0.5, sampletime=10,)
    @ntuple sima simb
end
safesave(simulationdir("sim2.jld2"), Dict("sim" => sim2))

## Simulation 3:
# as in simulation 2 but smaller proportion detected
sim3 = let 
    rng = Xoshiro(3)
    rnga = Xoshiro(3)
    rngb = Xoshiro(3)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:2]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:2]
    eta = 0.2
    sigma = 0.5
    phi = 0.05
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
    sima = packsimulations(rnga, 100, s1a, s2a; id="sim3a", minvalue=0.5, sampletime=10,)
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=50,
    )
    simb = packsimulations(rngb, 100, s1b, s2b; id="sim3b", minvalue=0.5, sampletime=10,)
    @ntuple sima simb
end
safesave(simulationdir("sim3.jld2"), Dict("sim" => sim3))

## Simulation 4:
# as in simulation 3 but with the intervention after the epidemic peak
sim4 = let 
    rng = Xoshiro(4)
    rnga = Xoshiro(4)
    rngb = Xoshiro(4)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:2]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:2]
    eta = 0.2
    sigma = 0.5
    phi = 0.05
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta2_a = x -> beta2ratio * beta1(x)
    beta2_b = x -> beta2_a(x) * (x < 80 ? 1 : 0.8)
    betas_a = [beta1, beta2_a]
    betas_b = [beta1, beta2_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=80,
    )
    sima = packsimulations(rnga, 130, s1a, s2a; id="sim4a", minvalue=0.5, sampletime=10,)
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=80,
    )
    simb = packsimulations(rngb, 130, s1b, s2b; id="sim4b", minvalue=0.5, sampletime=10,)
    @ntuple sima simb
end
safesave(simulationdir("sim4.jld2"), Dict("sim" => sim4))

## Simulation 5:
# as in simulation 3 but with four groups, two of which have the intervention
sim4 = let 
    rng = Xoshiro(5)
    rnga = Xoshiro(5)
    rngb = Xoshiro(5)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:4]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:4]
    eta = 0.2
    sigma = 0.5
    phi = 0.05
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta3ratio = rand(rng, Uniform(0.9, 1.1))
    beta4ratio = rand(rng, Uniform(0.9, 1.1))
    beta2 = x -> beta2ratio * beta1(x)
    beta3_a = x -> beta3ratio * beta1(x)
    beta3_b = x -> beta3_a(x) * (x < 40 ? 1 : 0.8)
    beta4_a = x -> beta4ratio * beta1(x)
    beta4_b = x -> beta4_a(x) * (x < 65 ? 1 : 0.8)
    betas_a = [beta1, beta2, beta3_a, beta4_a]
    betas_b = [beta1, beta2, beta3_b, beta4_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=nothing,
    )
    s3a = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_a[3], sigma, eta, phi, intervention=40,
    )
    s4a = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_a[4], sigma, eta, phi, intervention=65,
    )
    sima = packsimulations(
        rnga, 100, s1a, s2a, s3a, s4a; 
        id="sim5a", minvalue=0.5, sampletime=10,
    )
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=nothing,
    )
    s3b = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_b[3], sigma, eta, phi, intervention=40,
    )
    s4b = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_b[4], sigma, eta, phi, intervention=65,
    )
    simb = packsimulations(
        rngb, 100, s1b, s2b, s3b, s4b; 
        id="sim5b", minvalue=0.5, sampletime=10,
    )
    @ntuple sima simb
end
safesave(simulationdir("sim5.jld2"), Dict("sim" => sim5))

## Simulation 6:
# as in simulation 5 with competing interventions
sim6 = let 
    rng = Xoshiro(6)
    rnga = Xoshiro(6)
    rngb = Xoshiro(6)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:4]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:4]
    eta = 0.2
    sigma = 0.5
    phi = 0.05
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta3ratio = rand(rng, Uniform(0.9, 1.1))
    beta4ratio = rand(rng, Uniform(0.9, 1.1))
    beta2 = x -> beta2ratio * beta1(x) * (x < 75 ? 1 : 1.2)
    beta3_a = x -> beta3ratio * beta1(x)
    beta3_b = x -> beta3_a(x) * (x < 60 ? 1 : 0.8)
    beta4_a = x -> beta4ratio * beta1(x) * (x < 65 ? 1 : 1.2)
    beta4_b = x -> beta4_a(x) * (x < 45 ? 1 : 0.8)
    betas_a = [beta1, beta2, beta3_a, beta4_a]
    betas_b = [beta1, beta2, beta3_b, beta4_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=[nothing, 75],
    )
    s3a = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_a[3], sigma, eta, phi, intervention=[60, nothing],
    )
    s4a = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_a[4], sigma, eta, phi, intervention=[45, 65],
    )
    sima = packsimulations(
        rnga, 100, s1a, s2a, s3a, s4a; 
        id="sim6a", minvalue=0.5, sampletime=10,
    )
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=nothing,
    )
    s3b = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_b[3], sigma, eta, phi, intervention=40,
    )
    s4b = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_b[4], sigma, eta, phi, intervention=65,
    )
    simb = packsimulations(
        rngb, 100, s1b, s2b, s3b, s4b; 
        id="sim6b", minvalue=0.5, sampletime=10,
    )
    @ntuple sima simb
end
safesave(simulationdir("sim6.jld2"), Dict("sim" => sim6))

## Simulation 7:
# as in simulation 6 with gradually diverging transmission regardless of interventions
sim7 = let 
    rng = Xoshiro(7)
    rnga = Xoshiro(7)
    rngb = Xoshiro(7)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:4]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:4]
    eta = 0.2
    sigma = 0.5
    phi = 0.05
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta3ratio = rand(rng, Uniform(0.9, 1.1))
    beta4ratio = rand(rng, Uniform(0.9, 1.1))
    beta2 = x -> beta2ratio * beta1(x) * (x < 75 ? 1 : 1.2) * 0.002 * x
    beta3_a = x -> beta3ratio * beta1(x)
    beta3_b = x -> beta3_a(x) * (x < 60 ? 1 : 0.8) * 0.002 * x
    beta4_a = x -> beta4ratio * beta1(x) * (x < 65 ? 1 : 1.2)
    beta4_b = x -> beta4_a(x) * (x < 45 ? 1 : 0.8)
    betas_a = [beta1, beta2, beta3_a, beta4_a]
    betas_b = [beta1, beta2, beta3_b, beta4_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=[nothing, 75],
    )
    s3a = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_a[3], sigma, eta, phi, intervention=[60, nothing],
    )
    s4a = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_a[4], sigma, eta, phi, intervention=[45, 65],
    )
    sima = packsimulations(
        rnga, 100, s1a, s2a, s3a, s4a; 
        id="sim7a", minvalue=0.5, sampletime=10,
    )
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=nothing,
    )
    s3b = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_b[3], sigma, eta, phi, intervention=40,
    )
    s4b = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_b[4], sigma, eta, phi, intervention=65,
    )
    simb = packsimulations(
        rngb, 100, s1b, s2b, s3b, s4b; 
        id="sim7b", minvalue=0.5, sampletime=10,
    )
    @ntuple sima simb
end
safesave(simulationdir("sim7.jld2"), Dict("sim" => sim7))

## Simulation 8:
# as in simulation 6 with time-dependent proportion detected
sim8 = let 
    rng = Xoshiro(8)
    rnga = Xoshiro(8)
    rngb = Xoshiro(8)
    Ns = [round(Int, rand(rng, Uniform(500_000, 2_000_000))) for _ in 1:4]
    Es = rand.(rng, Binomial.(Ns, 0.0001))
    u0s = [simulationu0(S=(Ns[i] - Es[i]), E=Es[i]) for i in 1:4]
    eta = 0.2
    sigma = 0.5
    phi = x -> 0.04 + 0.0002 * x
    beta1 = x -> 0.4 + 0.02 * cos(x * 2pi / 365)
    beta2ratio = rand(rng, Uniform(0.9, 1.1))
    beta3ratio = rand(rng, Uniform(0.9, 1.1))
    beta4ratio = rand(rng, Uniform(0.9, 1.1))
    beta2 = x -> beta2ratio * beta1(x) * (x < 75 ? 1 : 1.2)
    beta3_a = x -> beta3ratio * beta1(x)
    beta3_b = x -> beta3_a(x) * (x < 60 ? 1 : 0.8)
    beta4_a = x -> beta4ratio * beta1(x) * (x < 65 ? 1 : 1.2)
    beta4_b = x -> beta4_a(x) * (x < 45 ? 1 : 0.8)
    betas_a = [beta1, beta2, beta3_a, beta4_a]
    betas_b = [beta1, beta2, beta3_b, beta4_b]
    s1a = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_a[1], sigma, eta, phi, intervention=nothing,
    )
    s2a = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_a[2], sigma, eta, phi, intervention=[nothing, 75],
    )
    s3a = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_a[3], sigma, eta, phi, intervention=[60, nothing],
    )
    s4a = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_a[4], sigma, eta, phi, intervention=[45, 65],
    )
    sima = packsimulations(
        rnga, 100, s1a, s2a, s3a, s4a; 
        id="sim8a", minvalue=0.5, sampletime=10,
    )
    s1b = packsimulationtuple( ; 
        u0=u0s[1], beta=betas_b[1], sigma, eta, phi, intervention=nothing,
    )
    s2b = packsimulationtuple( ; 
        u0=u0s[2], beta=betas_b[2], sigma, eta, phi, intervention=nothing,
    )
    s3b = packsimulationtuple( ; 
        u0=u0s[3], beta=betas_b[3], sigma, eta, phi, intervention=40,
    )
    s4b = packsimulationtuple( ; 
        u0=u0s[4], beta=betas_b[4], sigma, eta, phi, intervention=65,
    )
    simb = packsimulations(
        rngb, 100, s1b, s2b, s3b, s4b; 
        id="sim8b", minvalue=0.5, sampletime=10,
    )
    @ntuple sima simb
end
safesave(simulationdir("sim8.jld2"), Dict("sim" => sim8))
