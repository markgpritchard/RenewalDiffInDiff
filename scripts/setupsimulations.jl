
# To aid reproducibility, simulations are stored in the folder `data/sims`. This code was
# used to generate those simulations but it is not anticipated that it will need to be used
# again.

using DrWatson 
@quickactivate :RenewalDiffInDiff
include(srcdir("AnalysisFunctions.jl"))
using .AnalysisFunctions
include("simulationtransmissionparameters.jl")
using Random 

Random.seed!(1729)

# Two locations and two discrete transmission parameters 
simulation1dataset = let  
    interventions = InterventionsMatrix([ nothing, 100 ], 200)
    Ns = [ 8_000_000, 5_000_000 ]
    betafunctions = [ betafunction1_1, betafunction1_2 ]
    betafunctions_counterfactual = [ betafunction1_1, betafunction1_2_counterfactual ]
    seirparameters1(β) = SEIRParameters(β, 1 / 2, 1 / 2.5, 0.8)

    u0_1_1 = SEIRCompartments(Ns[1] - 500, 500)
    p_1_1 = seirparameters1(betafunction1_1)
    Random.seed!(11)
    sim_1_1 = runseir_noisy(u0_1_1, p_1_1, 200).reportedcases
    
    u0_1_2 = SEIRCompartments(Ns[2] - 50, 50)
    p_1_2 = seirparameters1(betafunction1_2)
    Random.seed!(12)
    sim_1_2 = runseir_noisy(u0_1_2, p_1_2, 200).reportedcases

    p_1_2_counterfactual = seirparameters1(betafunction1_2_counterfactual)
    Random.seed!(12)
    sim_1_2_counterfactual = runseir_noisy(u0_1_2, p_1_2_counterfactual, 200).reportedcases

    cases = hcat(sim_1_1, sim_1_2)
    cases_counterfactual = hcat(sim_1_1, sim_1_2_counterfactual)

    @ntuple betafunctions betafunctions_counterfactual cases cases_counterfactual interventions Ns
end

safesave(datadir("sims", "simulation1dataset.jld2"), ntuple2dict(simulation1dataset))

simulation2dataset = let  
    interventions = InterventionsMatrix([ nothing, 75, 100 ], 200)
    Ns = [ 18_000_000, 4_000_000, 20_000_000 ]
    betafunctions = [ betafunction2_1, betafunction2_2, betafunction2_3 ]
    betafunctions_counterfactual = [ 
        betafunction2_1, betafunction2_2_counterfactual, betafunction2_3_counterfactual 
    ]
    seirparameters2(β) = SEIRParameters(β, 1 / 2, 1 / 2.5, 0.45)

    u0_2_1 = SEIRCompartments(Ns[1] - 500, 500)
    p_2_1 = seirparameters2(betafunction2_1)
    Random.seed!(21)
    sim_2_1 = runseir_noisy(u0_2_1, p_2_1, 200).reportedcases
    
    u0_2_2 = SEIRCompartments(Ns[2] - 50, 50)
    p_2_2 = seirparameters2(betafunction2_2)
    Random.seed!(22)
    sim_2_2 = runseir_noisy(u0_2_2, p_2_2, 200).reportedcases

    u0_2_3 = SEIRCompartments(Ns[3] - 10, 10)
    p_2_3 = seirparameters2(betafunction2_3)
    Random.seed!(23)
    sim_2_3 = runseir_noisy(u0_2_3, p_2_3, 200).reportedcases

    p_2_2_counterfactual = seirparameters2(betafunction2_2_counterfactual)
    Random.seed!(22)
    sim_2_2_counterfactual = runseir_noisy(u0_2_2, p_2_2_counterfactual, 200).reportedcases

    p_2_3_counterfactual = seirparameters2(betafunction2_3_counterfactual)
    Random.seed!(23)
    sim_2_3_counterfactual = runseir_noisy(u0_2_3, p_2_3_counterfactual, 200).reportedcases

    cases = hcat(sim_2_1, sim_2_2, sim_2_3)
    cases_counterfactual = hcat(sim_2_1, sim_2_2_counterfactual, sim_2_3_counterfactual)

    @ntuple betafunctions betafunctions_counterfactual cases cases_counterfactual interventions Ns
end

safesave(datadir("sims", "simulation2dataset.jld2"), ntuple2dict(simulation2dataset))


simulation3dataset = let  
    interventions = InterventionsMatrix([ nothing, nothing, 75, 125 ], 200)
    secondaryinterventions = InterventionsMatrix([ nothing, 110, nothing, 50 ], 200)
    Ns = [ 8_000_000, 24_000_000, 13_000_000, 16_000_000 ]

    betafunctions = [ betafunction3_1, betafunction3_2, betafunction3_3, betafunction3_4 ]
    betafunctions_counterfactual = [ 
        betafunction3_1, betafunction3_2, 
        betafunction3_3_counterfactual, betafunction3_4_counterfactual
    ]
    seirparameters3(β) = SEIRParameters(β, 1 / 2, 1 / 2.5, 0.3)

    u0_3_1 = SEIRCompartments(Ns[1] - 5000, 5000)
    p_3_1 = seirparameters3(betafunction3_1)
    Random.seed!(31)
    sim_3_1 = runseir_noisy(u0_3_1, p_3_1, 200).reportedcases
    
    u0_3_2 = SEIRCompartments(Ns[2] - 8000, 8000)
    p_3_2 = seirparameters3(betafunction3_2)
    Random.seed!(32)
    sim_3_2 = runseir_noisy(u0_3_2, p_3_2, 200).reportedcases

    u0_3_3 = SEIRCompartments(Ns[3] - 1000, 1000)
    p_3_3 = seirparameters3(betafunction3_3)
    Random.seed!(33)
    sim_3_3 = runseir_noisy(u0_3_3, p_3_3, 200).reportedcases

    u0_3_4 = SEIRCompartments(Ns[4] - 10_000, 10_000)
    p_3_4 = seirparameters3(betafunction3_4)
    Random.seed!(34)
    sim_3_4 = runseir_noisy(u0_3_4, p_3_4, 200).reportedcases

    p_3_3_counterfactual = seirparameters3(betafunction3_3_counterfactual)
    Random.seed!(33)
    sim_3_3_counterfactual = runseir_noisy(u0_3_3, p_3_3_counterfactual, 200).reportedcases

    p_3_4_counterfactual = seirparameters3(betafunction3_4_counterfactual)
    Random.seed!(34)
    sim_3_4_counterfactual = runseir_noisy(u0_3_4, p_3_4_counterfactual, 200).reportedcases

    cases = hcat(sim_3_1, sim_3_2, sim_3_3, sim_3_4)
    cases_counterfactual = hcat(
        sim_3_1, sim_3_2, sim_3_3_counterfactual, sim_3_4_counterfactual
    )

    @ntuple betafunctions betafunctions_counterfactual cases cases_counterfactual interventions Ns secondaryinterventions
end

safesave(datadir("sims", "simulation3dataset.jld2"), ntuple2dict(simulation3dataset))

simulation4dataset = let  
    interventions = InterventionsMatrix([ nothing, 75, 100 ], 200)
    Ns = [ 18_000_000, 4_000_000, 20_000_000 ]

    betafunctions = [ betafunction4_1, betafunction4_2, betafunction4_3 ]
    betafunctions_counterfactual = [ 
        betafunction4_1, betafunction4_2_counterfactual, betafunction4_3_counterfactual 
    ]
    seirparameters4(β) = SEIRParameters(β, 1 / 2, 1 / 2.5, 0.45)

    u0_4_1 = SEIRCompartments(Ns[1] - 500, 500)
    p_4_1 = seirparameters4(betafunction4_1)
    Random.seed!(41)
    sim_4_1 = runseir_noisy(u0_4_1, p_4_1, 200).reportedcases
    
    u0_4_2 = SEIRCompartments(Ns[2] - 50, 50)
    p_4_2 = seirparameters4(betafunction4_2)
    Random.seed!(42)
    sim_4_2 = runseir_noisy(u0_4_2, p_4_2, 200).reportedcases

    u0_4_3 = SEIRCompartments(Ns[3] - 10, 10)
    p_4_3 = seirparameters4(betafunction4_3)
    Random.seed!(43)
    sim_4_3 = runseir_noisy(u0_4_3, p_4_3, 200).reportedcases

    p_4_2_counterfactual = seirparameters4(betafunction4_2_counterfactual)
    Random.seed!(42)
    sim_4_2_counterfactual = runseir_noisy(u0_4_2, p_4_2_counterfactual, 200).reportedcases

    p_4_3_counterfactual = seirparameters4(betafunction4_3_counterfactual)
    Random.seed!(43)
    sim_4_3_counterfactual = runseir_noisy(u0_4_3, p_4_3_counterfactual, 200).reportedcases

    cases = hcat(sim_4_1, sim_4_2, sim_4_3)
    cases_counterfactual = hcat(sim_4_1, sim_4_2_counterfactual, sim_4_3_counterfactual)

    @ntuple betafunctions betafunctions_counterfactual cases cases_counterfactual interventions Ns
end

safesave(datadir("sims", "simulation4dataset.jld2"), ntuple2dict(simulation4dataset))
