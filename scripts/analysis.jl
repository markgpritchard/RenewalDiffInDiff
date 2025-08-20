
using DrWatson
@quickactivate :RenewalDiffInDiff

using Distributions
using RenewalDiD

id = ARGS[1]

sim = load(simulationdir("sim$id.jld2"))["sim"]
model1 = renewaldid(                      
    sim1, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2), 1), 
        mu_delayprior=log(5),
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(8, 2),
        tauprior=Normal(0, 0.2),
    );                          
    mu=0.2, kappa=0.5,               
)

test = analysisworkflow(
    model1; 
    name="test1", 
    chain=1, 
    npriors=1000, 
    mapmaxtime=60, 
    nsamples=12, 
    priorsseed=1, 
    sampleseed=1001,
)

#test 

#priorsdf, priorschain = priorsworkflow(model1; chain=1, name="test1", npriors=1000, priorsseed=1)
#map_df, map_estimate = maximumlikelihoodworkflow(model1, priorsdf; chain=1, name="test1")


#=
test = analysisworkflow(
    model1; 
    name="test1", 
    chain=1, 
    npriors=1000, 
    mapmaxtime=60, 
    nsamples=12, 
    priorsseed=1, 
    sampleseed=1001,
)
=#

test = load(datadir("sims", "test1_map_1.jld2"))

#priorschain1 = sample(rng, model1, Prior(), 10_000)
#priorsdf1 = DataFrame(priorschain1)
priorsdf1 = test["priorsdf"]
priortraceplot1 = trplot(priorsdf1; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
priorsfittedoutputs1 = samplerenewaldidinfections(
    g_seir, priorsdf1, sim1; 
    mu=0.2, kappa=0.5,   
)
priorsoutputquantiles1 = quantilerenewaldidinfections(
    priorsfittedoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot1 = plotmodel(
    priorsoutputquantiles1, sim1; 
    interventionlinestyle=(:dot, :dense), linewidth=1,
)

#=
initindices1 = findall(x -> x <= 4, ordinalrank(priorsdf1.lp; rev=true)) 
priorsfittedinitoutputs1 = samplerenewaldidinfections(
    g_seir, priorsdf1, sim1, initindices1; 
    gamma=0.2, sigma=0.5,
)
priorsoutputinitquantiles1 = quantilerenewaldidinfections(
    priorsfittedinitoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsinitplot1 = plotmodel(
    priorsoutputinitquantiles1, sim1; 
    interventionlinestyle=(:dot, :dense), linewidth=1,
)
=#
#priorinitparams = [[values(priorsdf[i, 3:735])...] for i in initindices]
#priorinitparams1 = [[values(priorsdf1[i, 3:733])...] for i in initindices1]

mapdf = test["mapdf"]
mapoutputs1 = samplerenewaldidinfections(
    g_seir, mapdf, sim1; 
    mu=0.2, kappa=0.5, repeatsamples=1000,  
)
mapoutputquantiles1 = quantilerenewaldidinfections(
    mapoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
mapplot1 = plotmodel(
    mapoutputquantiles1, sim1; 
    interventionlinestyle=(:dot, :dense), linewidth=1,
)

#chain1 = sample(
#    rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 100, 4; 
#    initial_params=priorinitparams1,
#) 
#df1 = DataFrame(chain1)
df1 = test["mcmcdf"]
p1 = trplot(df1; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
#p2 = tracerankplot(df1; binsize=5, ncols=5, nplots=50, size=(1000, 1000)) 

fittedoutputs1 = samplerenewaldidinfections(
    g_seir, df1, sim1; 
    mu=0.2, kappa=0.5,   
)
outputquantiles1 = quantilerenewaldidinfections(
    fittedoutputs1, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
p3 = plotmodel(
    outputquantiles1, sim1;
    interventionlinestyle=(:dot, :dense), linewidth=1,
)



