
using DrWatson 
@quickactivate :RenewalDiffInDiff
using Random, StochasticTransitionModels 

function seirrates(u, t, p)
    s, e, i, i′, r = u  # i′ represents diagnosed infections. i + i′ is the total infectiouse prevalence
    n = sum(@view u[1:5])  # 6th compartment is cumulative diagnosed infecitons
    return [
        p.β(t) * s * (i + i′) / n,  # infection rate
        p.μ * e,  # end of latent period 
        p.θ * i,  # diagnosis 
        p.γ * i,  # recovery (undiagnosed)
        p.γ * i′  # recovery (diagnosed)
    ]
end

seirtransitionmatrix = [
    # s   e   i   i′  r   cumulative 
     -1   1   0   0   0   0    # infection rate
      0  -1   1   0   0   0    # end of latent period 
      0   0  -1   1   0   1    # diagnosis 
      0   0  -1   0   1   0    # recovery (undiagnosed)
      0   0   0  -1   1   0    # recovery (diagnosed)
]

Random.seed!(1729)

# Two locations and two discrete transmission parameters 

sim1parameters(beta) = SEIRParameters(beta, 0.5, 0.4, 0.8)
beta1a(t) = t <= 50 ? 0.6 : 0.66
beta1bcounterfactual(t) = 1.15 * beta1a(t)
beta1b(t) = t <= 50 ? beta1bcounterfactual(t) : 0.8 * beta1bcounterfactual(t)

simulation1dataset = let  
    interventions = InterventionsMatrix([ nothing, 50 ], 100)
    
    u01a = [ 8_000_000 - 750, 750, 0, 0, 0, 0  ]
    p1a = sim1parameters(beta1a)
    sim1a = stochasticmodel(seirrates, u01a, 1:100, p1a, seirtransitionmatrix)
    
    u01b = [ 5_000_000 - 200, 200, 0, 0, 0, 0 ]
    p1bcounterfactual = sim1parameters(beta1bcounterfactual)
    sim1bcounterfactual = stochasticmodel(
        seirrates, u01b, 1:100, p1bcounterfactual, seirtransitionmatrix
    )
    
    p1b = sim1parameters(beta1b)
    sim1b = vcat(
        sim1bcounterfactual[1:49, :],
        stochasticmodel(
            seirrates, sim1bcounterfactual[50, :], 50:100, p1b, seirtransitionmatrix
        )
    )
    
    prevalence = hcat(sim1a[:, 4], sim1b[:, 4])
    counterfactualprevalence = hcat(sim1a[:, 4], sim1bcounterfactual[:, 4])

    cases = zeros(Int, 100, 2)
    for t ∈ 2:100 
        cases[t, 1] = sim1a[t, 6] - sim1a[t-1, 6]
        cases[t, 2] = sim1b[t, 6] - sim1b[t-1, 6]
    end

    counterfactualcases = zeros(Int, 100, 2)
    for t ∈ 2:100 
        counterfactualcases[t, 1] = sim1a[t, 6] - sim1a[t-1, 6]
        counterfactualcases[t, 2] = sim1bcounterfactual[t, 6] - sim1bcounterfactual[t-1, 6]
    end

    @ntuple cases counterfactualcases interventions prevalence counterfactualprevalence Ns=[ 8_000_000, 5_000_000 ]
end

# Three locations, continuously changing transmission parameters, and a competing intervention  

sim2parameters(beta) = SEIRParameters(beta, 0.5, 0.4, 0.5)
beta2a(t) = 0.5 + 0.15 * cos(2π * (t - 80) / 365)
beta2bcounterfactual(t) = 1.15 * beta2a(t)
beta2b(t) = t <= 50 ? beta2bcounterfactual(t) : 0.8 * beta2bcounterfactual(t)
beta2ccounterfactual(t) = t <= 30 ? 0.9 * beta2a(t) : 0.9 * 1.15 * beta2a(t)
beta2c(t) = t <= 70 ? beta2ccounterfactual(t) : 0.8 * beta2ccounterfactual(t)

simulation2dataset = let  
    interventions = InterventionsMatrix([ nothing, 50, 70 ], 100)
    
    u02a = [ 7_000_000 - 1000, 1000, 0, 0, 0, 0  ]
    p2a = sim2parameters(beta2a)
    sim2a = stochasticmodel(seirrates, u02a, 1:100, p2a, seirtransitionmatrix)

    u02b = [ 6_000_000 - 2000, 2000, 0, 0, 0, 0 ]
    p2bcounterfactual = sim2parameters(beta2bcounterfactual)
    sim2bcounterfactual = stochasticmodel(
        seirrates, u02b, 1:100, p2bcounterfactual, seirtransitionmatrix
    )

    p2b = sim2parameters(beta2b)
    sim2b = vcat(
        sim2bcounterfactual[1:49, :],
        stochasticmodel(
            seirrates, sim2bcounterfactual[50, :], 50:100, p2b, seirtransitionmatrix
        )
    )

    u02c = [ 8_000_000 - 2000, 2000, 0, 0, 0, 0 ]
    p2ccounterfactual = sim2parameters(beta2ccounterfactual)
    sim2ccounterfactual = stochasticmodel(
        seirrates, u02c, 1:100, p2ccounterfactual, seirtransitionmatrix
    )

    p2c = sim2parameters(beta2c)
    sim2c = vcat(
        sim2ccounterfactual[1:69, :],
        stochasticmodel(
            seirrates, sim2ccounterfactual[70, :], 70:100, p2c, seirtransitionmatrix
        )
    )

    prevalence = hcat(sim2a[:, 4], sim2b[:, 4], sim2c[:, 4])
    counterfactualprevalence = hcat(
        sim2a[:, 4], sim2bcounterfactual[:, 4], sim2ccounterfactual[:, 4]
    )

    cases = zeros(Int, 100, 3)
    for t ∈ 2:100 
        cases[t, 1] = sim2a[t, 6] - sim2a[t-1, 6]
        cases[t, 2] = sim2b[t, 6] - sim2b[t-1, 6]
        cases[t, 3] = sim2c[t, 6] - sim2c[t-1, 6]
    end

    counterfactualcases = zeros(Int, 100, 3)
    for t ∈ 2:100 
        cases[t, 1] = sim2a[t, 6] - sim2a[t-1, 6]
        cases[t, 2] = sim2bcounterfactual[t, 6] - sim2bcounterfactual[t-1, 6]
        cases[t, 3] = sim2ccounterfactual[t, 6] - sim2ccounterfactual[t-1, 6]
    end

    @ntuple cases counterfactualcases interventions prevalence counterfactualprevalence Ns=[ 7_000_000, 6_000_000, 8_000_000 ]
end






