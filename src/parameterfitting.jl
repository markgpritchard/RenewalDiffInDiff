
function generatew_gt(f, cases, Ns; blankn=0)
    testf(f)
    w_gt = Matrix{Float64}(undef, size(cases))
    for g ∈ axes(w_gt, 2)
        for t ∈ axes(w_gt, 1)
            if t <= blankn 
                w_gt[t, g] = NaN
            elseif t == 1 
                w_gt[t, g] = _generatew_gtrow_cases(cases, t, g)
            else
                w_gt[t, g] = +(
                    _generatew_gtrow_cases(cases, t, g),
                    _generatew_gtrow_gentime(f, cases, t, g)
                ) 
            end
        end
    end
    return w_gt
end

_generatew_gtrow_cases(cases, t, g) = log(cases[t, g])

function _generatew_gtrow_gentime(f::Function, cases, t, g)
    return -log(sum([ f(τ) * cases[t-τ, g] for τ ∈ 1:(t - 1) ]))
end

function _generatew_gtrow_gentime(f::AbstractVector, cases, t, g)
    taumax = min(t - 1, length(f))
    return -log(sum([ f[τ] * cases[t-τ, g] for τ ∈ 1:taumax ]))
end

@model _estimatelogzeta_g(μ, σ2) = logzeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogeta_t(μ, σ2) = logeta ~ Normal(μ, sqrt(σ2))
@model _estimatelogeta_t(etaprior) = logeta ~ etaprior
@model _estimatelogsecondarydelta(deltaprior) = logsecondarydelta ~ deltaprior
_elsd(deltaprior) = _estimatelogsecondarydelta(deltaprior)

function renewaldiffindiffparameters(; w, y, interventions, timeknots, Ns, kwargs...)
    return renewaldiffindiffparameters(w, y, interventions, timeknots, Ns; kwargs...)
end

@model function renewaldiffindiffparameters(
    w, y, interventions, timeknots, Ns;
    etavarianceprior=Exponential(1),
    extrapl=zeros(1), extrapr=zeros(1), 
    interceptprior=Normal(0, 1),
    logdeltaprior=Normal(0, 1),
    psiprior=automatic,
    secondaryinterventions=nothing,
    sigmasquareprior=Exponential(1),
    zeroperiod=2,
    zetavarianceprior=Exponential(1),
)
    ntimeknots = length(timeknots)
    nlocations = size(w, 2)
    psiminimum = maximum([ sum(y[:, g]) / Ns[g] for g ∈ axes(w, 2) ])
    @assert 0 <= psiminimum <= 1

    eta_sigma2 ~ etavarianceprior
    intercept ~ interceptprior
    logdelta ~ logdeltaprior
    if typeof(psiprior) == Automatic 
        psi ~ Uniform(psiminimum, 1) 
    else
        psi ~ truncated(psiprior, psiminimum, 1) 
    end
    σ2 ~ sigmasquareprior 
    zeta_sigma2 ~ zetavarianceprior

    logzeta_g_vec = Vector{typeof(zeta_sigma2)}(undef, nlocations) 
    for i ∈ 1:nlocations
        @submodel prefix="logzeta_g$i" logzeta = _estimatelogzeta_g(0, 1)
        logzeta_g_vec[i] = intercept + logzeta * zeta_sigma2
    end

    logeta_t_vec = Vector{typeof(eta_sigma2)}(undef, ntimeknots) 
    for i ∈ 1:ntimeknots
        if i == zeroperiod
            logeta = zero(typeof(eta_sigma2))
        else
            @submodel prefix="eta_t$i" logeta = _estimatelogeta_t(0, 1)
        end
        logeta_t_vec[i] = logeta * eta_sigma2
    end
    timespline = CubicSpline(timeknots, ForwardDiff.value.(logeta_t_vec); extrapl, extrapr)

    if isnothing(secondaryinterventions)
        logsecondarydelta = nothing 
    elseif isa(secondaryinterventions, AbstractMatrix)
        logsecondarydelta ~ logdeltaprior
    else 
        logsecondarydelta = Vector{typeof(logdelta)}(undef, length(secondaryinterventions)) 
        for i ∈ eachindex(secondaryinterventions)
            @submodel prefix="logsecondarydelta$i" logsecondarydeltai = _elsd(logdeltaprior)
            logsecondarydelta[i] = logsecondarydeltai
        end
    end

    _Mopt = typeof(
        logeffectivereproductionratio(
            logzeta_g_vec,
            timespline, 
            logdelta,
            interventions, 
            0,
            logsecondarydelta,
            secondaryinterventions,
            1, 
            1
        )
    )
    modeloutput = Matrix{_Mopt}(undef, size(w))
    Threads.@threads for g ∈ axes(w, 2)
        for t ∈ axes(w, 1)
            s = _calcs(y, Ns, psi, t, g)
            modeloutput[t, g] = logeffectivereproductionratio(
                logzeta_g_vec[g],
                timespline,
                logdelta,
                interventions, 
                s,
                logsecondarydelta,
                secondaryinterventions,
                t, 
                g
            )
        end
    end

    for i ∈ eachindex(w)
        _skip(w[i]) && continue 
        w[i] ~ Normal(modeloutput[i], σ2)
    end 
end

_skip(x) = isnan(x) || x == -Inf || x == Inf 

function loadanalysisdictsasdf(modelname, n_chains, maxrounds, seedstart; nsims=4)
    fn1 = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + 1)
    chain1 = load(fn1)["chain"]
    df = DataFrame(chain1)
    for i ∈ 2:nsims 
        fn = _findanalysisfilename(modelname, n_chains, maxrounds, seedstart + i)
        isnothing(fn) && continue
        chain = load(fn)["chain"]
        _tdf = DataFrame(chain)
        _tdf.chain = [ i for _ ∈ axes(_tdf, 1) ]
        df = vcat(df, _tdf) 
    end
    return df
end

function _findanalysisfilename(modelname, n_chains, maxrounds, seed)
    n_rounds = maxrounds
    while n_rounds > 0 
        fn = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds)_seed=$(seed).jld2"
        if isfile(datadir("sims", fn))
            return datadir("sims", fn)
        else 
            n_rounds += -1 
        end
    end
    return nothing
end

function keyvalues(fitteddf, fittedvaluesset; cri=[ 0.05, 0.95 ])
    deltamean = exp(mean(fitteddf.logdelta))
    deltap05_95 = exp.(quantile(fitteddf.logdelta, cri))
    totalcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcases = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec),
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    peakcasesdate = Matrix{Int}(
        undef, 
        length(fittedvaluesset.y_matrix_poisson_vec), 
        size(fittedvaluesset.y_matrix_poisson_vec[1], 2)
    )
    for i ∈ eachindex(fittedvaluesset.y_matrix_poisson_vec)
        totalcases[i, :] = [
            sum(fittedvaluesset.y_matrix_poisson_vec[i][:, j]) - 
                sum(fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j])
            for j ∈ axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
        ]
        for j ∈ axes(fittedvaluesset.y_matrix_poisson_vec[i], 2)
            x, ind = findmax(fittedvaluesset.y_matrix_poisson_vec[i][:, j])
            x_cf, ind_cf = findmax(
                fittedvaluesset.y_matrix_poisson_vec_counterfactual[i][:, j]
            )
            peakcases[i, j] = x - x_cf 
            peakcasesdate[i, j] = ind - ind_cf
        end
    end
    totalcasesmeaneffect = [ mean(totalcases[:, i]) for i ∈ axes(totalcases, 2) ]
    totalcasesp05_90 = [ quantile(totalcases[:, i], cri) for i ∈ axes(totalcases, 2) ]
    peakcasesmeaneffect = [ mean(peakcases[:, i]) for i ∈ axes(peakcases, 2) ]
    peakcasesp05_90 = [ quantile(peakcases[:, i], cri) for i ∈ axes(peakcases, 2) ]
    peakcasesdatemeaneffect = [ mean(peakcasesdate[:, i]) for i ∈ axes(peakcasesdate, 2) ]
    peakcasesdatep05_90 = [ quantile(peakcasesdate[:, i], cri) for i ∈ axes(peakcasesdate, 2) ]
    return (
        deltamean=deltamean,
        deltap05_95=deltap05_95,
        totalcasesmeaneffect=totalcasesmeaneffect,
        totalcasesp05_90=totalcasesp05_90,
        peakcasesmeaneffect=peakcasesmeaneffect,
        peakcasesp05_90=peakcasesp05_90,
        peakcasesdatemeaneffect=peakcasesdatemeaneffect,
        peakcasesdatep05_90=peakcasesdatep05_90,
    )
end

