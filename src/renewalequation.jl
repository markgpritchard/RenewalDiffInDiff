
proportionwaned(t; mean=10) = pdf(Exponential(mean), t)

_notnegative(x::T) where T = max(x, zero(x))

function initialsusceptible(w, M, t, j)
    t >= 2 || throw(DomainError(t, "requires t ≥ 2"))
    return imag(M[t-1, j]) + sum([ real(M[x, j]) * w(t - x) for x ∈ 1:(t-1) ])
end

function expectedinfections(g, M, ρ, s, n, t, j; tmax=28) 
    t >= 2 || throw(DomainError(t, "requires t ≥ 2"))
    return s * ρ * sum([ real(M[x, j]) * g(t - x) for x ∈ max(1, t - tmax):t-1 ]) / n
end

@model function R0_did(
    incidence,
    interventions;
    g,
    N,
    seedtimes=7,
    alphaprior=Normal(log(2), 1),
    gammavarianceprior=Exponential(0.4),
    delaymeanprior=Exponential(1),
    delayvarianceprior=Exponential(1),
    firstthetaprior=Normal(0, 0.05),
    immunedurationprior=Exponential(100),
    rparameter_prior=Exponential(1),
    psiprior=Beta(6, 4),
    thetavarianceprior=Exponential(1), 
    tauprior=Normal(0, 1),
) 
    ntimes, nlocations = size(incidence)

    ## reproduction ratio 
    α ~ alphaprior

    gammavar ~ gammavarianceprior
    gammavec ~ filldist(Normal(0, 1), nlocations)

    θ1 ~ firstthetaprior
    thetavar ~ thetavarianceprior
    thetaothervec ~ filldist(Normal(0, 1), ntimes)
    thetavec = cumsum([ θ1; thetaothervec .* sqrt(thetavar) ])

    τ ~ tauprior

    Rmat = [
        exp(α + gammavec[j] * sqrt(gammavar) + thetavec[t] + interventions[t, j] * τ)
        for t ∈ 1:ntimes, j ∈ 1:nlocations
    ]

    ## infections
    meanimmuneduration ~ immunedurationprior
    w(x) = proportionwaned(x; mean=meanimmuneduration)

    Nmat ~ filldist(truncated(Normal(0, 1); lower=-1), ntimes + seedtimes, nlocations)

    Imat = Matrix{Complex{typeof(α)}}(undef, ntimes + seedtimes, nlocations)

    iota_zeros = [ 
        max(0, log(sum(@view incidence[1:seedtimes, j]) / seedtimes)) 
        for j ∈ 1:nlocations 
    ]
     
    for j ∈ 1:nlocations, t ∈ 1:seedtimes
        Imat[t, j] = (
            susceptible = t == 1 ? N[j] : initialsusceptible(w, Imat, t, j);
            exptdinfection = exp(iota_zeros[j] - 2 * (seedtimes - t));
            infection = min(susceptible, exptdinfection * (1 + Nmat[t, j]));
            complex(infection, susceptible - infection)
        )
    end

    for j ∈ 1:nlocations, t ∈ (seedtimes + 1):(seedtimes + ntimes)
        Imat[t, j] = (
            susceptible = _notnegative(initialsusceptible(w, Imat, t, j));
            exptdinfection = _notnegative(
                expectedinfections(
                    g, 
                    Imat, 
                    Rmat[t-seedtimes, j], 
                    susceptible, 
                    N[j],
                    t, 
                    j
                )
            );
            infection = min(susceptible, exptdinfection + sqrt(exptdinfection) * Nmat[t, j]);
            complex(_notnegative(infection), susceptible - _notnegative(infection)) 
        )
    end

    ## Delay in diagnosis 
    μ_D ~ delaymeanprior
    σ2_D ~ delayvarianceprior 
    θ_D = σ2_D / μ_D
    α_D = μ_D / θ_D

    Dmat = [
        sum([ Imat[x, j] * pdf(Gamma(α_D, θ_D), t - x) for x ∈ max(t-20, 1):(t-1) ])
        for t ∈ (seedtimes + 1):(seedtimes + ntimes), j ∈ 1:nlocations
    ]

    # proportion who are diagnosed 
    ψ ~ psiprior

    r ~ rparameter_prior
    #p = r ./ (r .+ real.(Imat[(seedtimes + 1):(seedtimes + ntimes), :]) .* ψ) 
    p = r ./ (r .+ real.(Dmat) .* ψ) 

    if r <= 0 || isnan(minimum(p)) || minimum(p) <= 0 || maximum(p) > 1 
        Turing.@addlogprob! -Inf
        return nothing
    end

    incidence ~ arraydist(NegativeBinomial.(r, p))
end
