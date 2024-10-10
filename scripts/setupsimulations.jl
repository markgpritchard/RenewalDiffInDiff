
using DrWatson 
@quickactivate :RenewalDiffInDiff
using Random, StochasticTransitionModels 

function seirrates(u, t, p::SEIRParameters{<:Function, <:Real, <:Real})
    s, e, i, i′, r = u  # i′ represents diagnosed infections. i + i′ is the total infectiouse prevalence
    n = sum(@view u[1:5])  # 6th compartment is cumulative diagnosed infecitons
    return [
        p.β(t) * s * (i + i′) / n,  # infection rate
        p.μ * e,  # end of latent period 
        p.θ * p.γ * i / (1 - p.θ),  # diagnosis 
        p.γ * i,  # recovery (undiagnosed)
        p.γ * i′  # recovery (diagnosed)
    ]
end

function seirrates(u, t, p::SEIRParameters{<:Function, <:Real, <:Function})
    s, e, i, i′, r = u  # i′ represents diagnosed infections. i + i′ is the total infectiouse prevalence
    n = sum(@view u[1:5])  # 6th compartment is cumulative diagnosed infecitons
    return [
        p.β(t) * s * (i + i′) / n,  # infection rate
        p.μ * e,  # end of latent period 
        p.θ(t) * p.γ * i / (1 - p.θ(t)),  # diagnosis 
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

simparameters(beta, theta) = SEIRParameters(beta, 0.5, 0.4, theta)

sim1parameters(beta) = simparameters(beta, 0.8)
beta1a(t) = t <= 50 ? 0.6 : 0.66
beta1bcounterfactual(t) = 1.15 * beta1a(t)
beta1b(t) = t <= 50 ? beta1bcounterfactual(t) : 0.8 * beta1bcounterfactual(t)

sim1a_parameters(beta) = simparameters(beta, 0.6)
beta1a_a(t) = 0.5 + 0.1 * cos(2π * (t - 20) / 365)
beta1a_bcounterfactual(t) = 1.15 * beta1a_a(t)
beta1a_b(t) = t <= 50 ? beta1a_bcounterfactual(t) : 0.8 * beta1a_bcounterfactual(t)

sim2parameters(beta) = simparameters(beta, 0.5)
beta2a(t) = 0.5 + 0.15 * cos(2π * (t - 80) / 365)
beta2bcounterfactual(t) = 1.15 * beta2a(t)
beta2b(t) = t <= 50 ? beta2bcounterfactual(t) : 0.8 * beta2bcounterfactual(t)
beta2ccounterfactual(t) = t <= 30 ? 0.9 * beta2a(t) : 0.9 * 1.15 * beta2a(t)
beta2c(t) = t <= 70 ? beta2ccounterfactual(t) : 0.8 * beta2ccounterfactual(t)

sim3parameters(beta) = simparameters(beta, 0.3)
beta3a(t) = 0.5 + 0.15 * cos(2π * (t - 80) / 365)
beta3bcounterfactual(t) = beta3a(t) * (0.9 + 0.005 * t)
beta3b(t) = t <= 50 ? beta3bcounterfactual(t) : 0.8 * beta3bcounterfactual(t)

sim4parameters(beta, theta) = simparameters(beta, theta)
beta4a(t) = beta2a(t)
beta4bcounterfactual(t) = beta2bcounterfactual(t)
beta4b(t) = beta2b(t)
theta4a(t) = 0.3
theta4b(t) = t <= 50 ? 0.3 : 1.2 * 0.3

if isfile(datadir("sims", "simulation1dataset.jld2"))
    simulation1dataset = load(datadir("sims", "simulation1dataset.jld2"))
else 
    simulation1dataset = let  
        interventions = InterventionsMatrix([ nothing, 50 ], 100)
        
        u01a = [ 400_000 - 40, 40, 0, 0, 0, 0  ]
        p1a = sim1parameters(beta1a)
        sim1a = stochasticmodel(seirrates, u01a, 1:100, p1a, seirtransitionmatrix)
        
        u01b = [ 250_000 - 10, 10, 0, 0, 0, 0 ]
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

        Dict(
            "cases" => cases, 
            "cases_counterfactual" => counterfactualcases,
            "interventions" => interventions, 
            "prevalence" => prevalence, 
            "counterfactualprevalence" => counterfactualprevalence, 
            "Ns" => [ 400_000, 250_000 ],
        )
    end

    safesave(datadir("sims", "simulation1dataset.jld2"), simulation1dataset)
end

# Two locations, continuously changing transmission parameters

if isfile(datadir("sims", "simulation1a_dataset.jld2"))
    simulation1a_dataset = load(datadir("sims", "simulation1a_dataset.jld2"))
else 
    simulation1a_dataset = let  
        interventions = InterventionsMatrix([ nothing, 50, 70 ], 100)
        
        u01a_a = [ 300_000 - 50, 50, 0, 0, 0, 0 ]
        p1a_a = sim1a_parameters(beta1a_a)
        sim1a_a = stochasticmodel(seirrates, u01a_a, 1:100, p1a_a, seirtransitionmatrix)

        u01a_b = [ 250_000 - 100, 100, 0, 0, 0, 0 ]
        p1a_bcounterfactual = sim1a_parameters(beta1a_bcounterfactual)
        sim1a_bcounterfactual = stochasticmodel(
            seirrates, u01a_b, 1:100, p1a_bcounterfactual, seirtransitionmatrix
        )

        p1a_b = sim1a_parameters(beta1a_b)
        sim1a_b = vcat(
            sim1a_bcounterfactual[1:49, :],
            stochasticmodel(
                seirrates, sim1a_bcounterfactual[50, :], 50:100, p1a_b, seirtransitionmatrix
            )
        )

        prevalence = hcat(sim1a_a[:, 4], sim1a_b[:, 4])
        counterfactualprevalence = hcat(sim1a_a[:, 4], sim1a_bcounterfactual[:, 4])

        cases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            cases[t, 1] = sim1a_a[t, 6] - sim1a_a[t-1, 6]
            cases[t, 2] = sim1a_b[t, 6] - sim1a_b[t-1, 6]
        end

        counterfactualcases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            counterfactualcases[t, 1] = sim1a_a[t, 6] - sim1a_a[t-1, 6]
            counterfactualcases[t, 2] = sim1a_bcounterfactual[t, 6] - sim1a_bcounterfactual[t-1, 6]
        end

        Dict(
            "cases" => cases, 
            "cases_counterfactual" => counterfactualcases,
            "interventions" => interventions, 
            "prevalence" => prevalence, 
            "counterfactualprevalence" => counterfactualprevalence, 
            "Ns" => [ 300_000, 250_000 ],
        )
    end

    safesave(datadir("sims", "simulation1a_dataset.jld2"), simulation1a_dataset)
end


# Three locations, continuously changing transmission parameters, and a competing intervention  

if isfile(datadir("sims", "simulation2dataset.jld2"))
    simulation2dataset = load(datadir("sims", "simulation2dataset.jld2"))
else 
    simulation2dataset = let  
        interventions = InterventionsMatrix([ nothing, 50, 70 ], 100)
        
        u02a = [ 350_000 - 50, 50, 0, 0, 0, 0 ]
        p2a = sim2parameters(beta2a)
        sim2a = stochasticmodel(seirrates, u02a, 1:100, p2a, seirtransitionmatrix)

        u02b = [ 300_000 - 100, 100, 0, 0, 0, 0 ]
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

        u02c = [ 400_000 - 100, 100, 0, 0, 0, 0 ]
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
            counterfactualcases[t, 1] = sim2a[t, 6] - sim2a[t-1, 6]
            counterfactualcases[t, 2] = sim2bcounterfactual[t, 6] - sim2bcounterfactual[t-1, 6]
            counterfactualcases[t, 3] = sim2ccounterfactual[t, 6] - sim2ccounterfactual[t-1, 6]
        end

        Dict(
            "cases" => cases, 
            "cases_counterfactual" => counterfactualcases,
            "interventions" => interventions, 
            "prevalence" => prevalence, 
            "counterfactualprevalence" => counterfactualprevalence, 
            "Ns" => [ 350_000, 300_000, 400_000 ],
        )
    end

    safesave(datadir("sims", "simulation2dataset.jld2"), simulation2dataset)
end

# Two locations, transmission parameters violate common trends  

if isfile(datadir("sims", "simulation3dataset.jld2"))
    simulation3dataset = load(datadir("sims", "simulation3dataset.jld2"))
else 
    simulation3dataset = let  
        interventions = InterventionsMatrix([ nothing, 50 ], 100)
        
        u03a = [ 250_000 - 100, 100, 0, 0, 0, 0 ]
        p3a = sim3parameters(beta3a)
        sim3a = stochasticmodel(seirrates, u03a, 1:100, p3a, seirtransitionmatrix)

        u03b = [ 250_000 - 100, 100, 0, 0, 0, 0 ]
        p3bcounterfactual = sim3parameters(beta3bcounterfactual)
        sim3bcounterfactual = stochasticmodel(
            seirrates, u03b, 1:100, p3bcounterfactual, seirtransitionmatrix
        )

        p3b = sim3parameters(beta3b)
        sim3b = vcat(
            sim3bcounterfactual[1:49, :],
            stochasticmodel(
                seirrates, sim3bcounterfactual[50, :], 50:100, p3b, seirtransitionmatrix
            )
        )

        prevalence = hcat(sim3a[:, 4], sim3b[:, 4])
        counterfactualprevalence = hcat(sim3a[:, 4], sim3bcounterfactual[:, 4])

        cases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            cases[t, 1] = sim3a[t, 6] - sim3a[t-1, 6]
            cases[t, 2] = sim3b[t, 6] - sim3b[t-1, 6]
        end

        counterfactualcases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            counterfactualcases[t, 1] = sim3a[t, 6] - sim3a[t-1, 6]
            counterfactualcases[t, 2] = sim3bcounterfactual[t, 6] - sim3bcounterfactual[t-1, 6]
        end

        Dict(
            "cases" => cases, 
            "cases_counterfactual" => counterfactualcases,
            "interventions" => interventions, 
            "prevalence" => prevalence, 
            "counterfactualprevalence" => counterfactualprevalence, 
            "Ns" => [ 250_000, 250_000 ],
        )
    end

    safesave(datadir("sims", "simulation3dataset.jld2"), simulation3dataset)
end

# Two locations, proportion detected changes at time of intervention

if isfile(datadir("sims", "simulation4dataset.jld2"))
    simulation4dataset = load(datadir("sims", "simulation4dataset.jld2"))
else 
    simulation4dataset = let  
        interventions = InterventionsMatrix([ nothing, 50 ], 100)
        
        u04a = [ 350_000 - 50, 50, 0, 0, 0, 0 ]
        p4a = sim4parameters(beta4a, theta4a)
        sim4a = stochasticmodel(seirrates, u04a, 1:100, p4a, seirtransitionmatrix)

        u04b = [ 150_000 - 50, 50, 0, 0, 0, 0 ]
        p4bcounterfactual = sim4parameters(beta4bcounterfactual, theta4b)
        sim4bcounterfactual = stochasticmodel(
            seirrates, u04b, 1:100, p4bcounterfactual, seirtransitionmatrix
        )

        p4b = sim4parameters(beta4b, theta4b)
        sim4b = vcat(
            sim4bcounterfactual[1:49, :],
            stochasticmodel(
                seirrates, sim4bcounterfactual[50, :], 50:100, p4b, seirtransitionmatrix
            )
        )

        prevalence = hcat(sim4a[:, 4], sim4b[:, 4])
        counterfactualprevalence = hcat(sim4a[:, 4], sim4bcounterfactual[:, 4])

        cases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            cases[t, 1] = sim4a[t, 6] - sim4a[t-1, 6]
            cases[t, 2] = sim4b[t, 6] - sim4b[t-1, 6]
        end

        counterfactualcases = zeros(Int, 100, 2)
        for t ∈ 2:100 
            counterfactualcases[t, 1] = sim4a[t, 6] - sim4a[t-1, 6]
            counterfactualcases[t, 2] = sim4bcounterfactual[t, 6] - sim4bcounterfactual[t-1, 6]
        end

        Dict(
            "cases" => cases, 
            "cases_counterfactual" => counterfactualcases,
            "interventions" => interventions, 
            "prevalence" => prevalence, 
            "counterfactualprevalence" => counterfactualprevalence, 
            "Ns" => [ 350_000, 150_000 ],
        )
    end

    safesave(datadir("sims", "simulation4dataset.jld2"), simulation4dataset)
end
