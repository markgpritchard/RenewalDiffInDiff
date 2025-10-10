
using DrWatson
@quickactivate :RenewalDiffInDiff

using Distributions
using RenewalDiD

id = parse(Int, ARGS[1])
chain = parse(Int, ARGS[2])
npriors = parse(Int, ARGS[3])
mapmaxtime = parse(Int, ARGS[4])
nsamples = parse(Int, ARGS[5])

#= use line below for arguments when running in REPL 
id = 1; chain = 1; npriors = 1000; mapmaxtime = 60; nsamples = 25;
=#

sampleseed = 1000 * id + chain

sim = let
    initsim = load(simulationdir("sim$id.jld2"))["sim"]
    offsetinterventions = addoffsetstointerventionarray(initsim.interventions)
    RenewalDiDData( ;
        observedcases=initsim.observedcases, 
        interventions=offsetinterventions, 
        Ns=initsim.Ns, 
        exptdseedcases=initsim.exptdseedcases, 
        id="$(initsim.id)leadlag"
    )
end

model = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(16, 4),
        tauprior=Normal(0, 0.2),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, 
    kappa=0.5,
    thetainterval=7,         
)

analysis = analysisworkflow(
    model; 
    name="analysis$(id)_offset", 
    chain, 
    npriors, 
    mapmaxtime, 
    nsamples, 
    priorsseed=id, 
    sampleseed,
)
