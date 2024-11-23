
using CairoMakie, DataFrames, Pigeons, PlotFormatting, Turing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", "chain", "lp", "n_steps", "is_accept", "acceptance_rate", "log_density", 
    "hamiltonian_energy", "hamiltonian_energy_error", "max_hamiltonian_energy_error", 
    "tree_depth", "numerical_error", "step_size", "nom_step_size"
]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting functions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotchains(data::Chains; kwargs...) = plotchains(DataFrame(data); kwargs...)

function plotchains(data::DataFrame; size=( 400, 5000 ), kwargs...)
    fig = Figure(; size)
    plotchains!(fig, data; kwargs...)
    return fig
end

function plotchains!(fig::Figure, data::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plotchains!(gl, data; kwargs...)
end

function plotchains!(
    gl::GridLayout, data::DataFrame; 
    colnames=RenewalDiffInDiff.automatic, plotnames_ind=RenewalDiffInDiff.automatic,
    ylabels=RenewalDiffInDiff.automatic, kwargs...
)
    @unpack colnames, plotnames_ind = processplotchains(
        data;
        colnames, plotnames_ind, kwargs...
    )
        
    ax = [ Axis(gl[i, 1]) for i ∈ eachindex(plotnames_ind) ]
    n_ax = length(ax)
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j], linewidth=1,)
            if j == 1 
                if ylabels == RenewalDiffInDiff.automatic
                    Label(gl[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
                else 
                    Label(gl[i, 0], ylabels[i]; tellheight=false)
                end
            end
        end
    end

    for (i, a) ∈ enumerate(ax)
        if i == n_ax 
            formataxis!(a)
        else 
            formataxis!(a; hidex=true)
        end
    end

    Label(gl[n_ax+1, 1], "sample"; fontsize=11.84, tellwidth=false)

    colgap!(gl, 1, 5)
    rowgap!(gl, n_ax, 5)
end

function _processplotchains(
    data, ::RenewalDiffInDiff.Automatic, ::RenewalDiffInDiff.Automatic; 
    logdensity="log_density"
)
    return _processplotchains(data; logdensity)
end

function _processplotchains(
    data, cn, ::RenewalDiffInDiff.Automatic; 
    logdensity="log_density"
)
    @unpack plotnames_ind = _processplotchains(data; kwargs...)
    return @ntuple colnames=cn plotnames_ind
end

function _processplotchains(
    data, ::RenewalDiffInDiff.Automatic, pni; 
    logdensity="log_density"
)
    @unpack colnames = _processplotchains(data; kwargs...)
    return @ntuple colnames plotnames_ind=pni

end

function _processplotchains(data, colnames, plotnames_ind; logdensity="log_density")
    return @ntuple colnames plotnames_ind
end

function _processplotchains(data; logdensity="log_density")
    # "log_density" is the label given by `Pigeons` output. Turing labels it "lp".
    colnames = names(data)
    lp_ind = findall(x -> x == logdensity, colnames)
    _plotnames_ind = findall(x -> x ∉ _NOPLOTNAMES, colnames)
    plotnames_ind = [ lp_ind; _plotnames_ind ]
    return @ntuple colnames plotnames_ind
end

function processplotchains(data; colnames=RenewalDiffInDiff.automatic, plotnames_ind=RenewalDiffInDiff.automatic, logdensity="log_density")
    _processplotchains(data, colnames, plotnames_ind; logdensity)
end

function plotrenewalequationsamples(args...; plotsize=( 800, 800 ), kwargs...)
    fig = Figure(; size=plotsize)
    plotrenewalequationsamples!(fig, args...; kwargs...)
    return fig
end

function plotrenewalequationsamples!(fig::Figure, args...; kwargs...)
    gl = GridLayout(fig[1, 1])
    plotrenewalequationsamples!(gl, args...; kwargs...)
end

function plotrenewalequationsamples!(
    gl::GridLayout, dataset::Dict, w, fittedvaluesset; 
    kwargs...
)
    @unpack cases, cases_counterfactual, Ns = dataset 
    fittedws = fitws(cases, Ns, fittedvaluesset)
    plotrenewalequationsamples!(
        gl, cases, cases_counterfactual, w, Ns, fittedvaluesset, fittedws; 
        kwargs...
    )
end

function plotrenewalequationsamples!(
    gl::GridLayout, dataset::NamedTuple, w, fittedvaluesset; 
    kwargs...
) 
    @unpack cases, counterfactualcases, Ns = dataset 
    fittedws = fitws(cases, Ns, fittedvaluesset)
    plotrenewalequationsamples!(
        gl, cases, counterfactualcases, w, Ns, fittedvaluesset, fittedws; 
        kwargs...
    )
end

function plotrenewalequationsamples!(
    gl::GridLayout, cases::Matrix, w, Ns::Vector, fittedvaluesset; 
    kwargs...
)
    fittedws = fitws(cases, Ns, fittedvaluesset)
    plotrenewalequationsamples!(gl, cases, w, Ns, fittedvaluesset, fittedws; kwargs...)
end

function plotrenewalequationsamples!(
    gl::GridLayout, cases::Matrix, w, Ns::Vector, fittedvaluesset, fittedws; 
    kwargs...
)
    plotrenewalequationsamples!(
        gl, cases, nothing, w, Ns, fittedvaluesset, fittedws; 
        kwargs...
    )
end

function plotrenewalequationsamples!(
    gl::GridLayout, 
    cases::Matrix, cases_counterfactual, w, Ns::Vector, fittedvaluesset, fittedws;
    betafunctions=nothing, betafunctions_counterfactual=nothing,
    datacolour=:black, simcolour=COLOURVECTOR[2], fittedcolour=( COLOURVECTOR[1], 0.75 ), 
    counterfactualcolour=( COLOURVECTOR[2], 0.75 ),
    infectiousduration=1, markersize=3, rhoclip=Inf,
    columntitles=nothing, columntitlefontsize=11.84, 
    plotrhocounterfactuals=false,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time"
)
    nlocations = size(cases, 2)
    for i ∈ 1:nlocations
        if !isnothing(columntitles)
            Label(
            gl[0, i], columntitles[i]; 
            fontsize=columntitlefontsize, halign=:left, tellwidth=false
        )
        end
    end

    plotrenewalequationsamples_w!(
        gl, cases, w, fittedvaluesset, fittedws, 1;
        datacolour, fittedcolour, markersize, hidex=true, xticklabelrotation, xticks, xtitle
    )
    plotrenewalequationsamples_r0!(
        gl, cases, fittedvaluesset, 2;
        betafunctions, simcolour, fittedcolour, infectiousduration, markersize, rhoclip,
        plotcounterfactuals=plotrhocounterfactuals, counterfactualcolour, 
        hidex=true, xticklabelrotation, xticks, xtitle
    )
    plotrenewalequationsamples_cases!(
        gl, cases, Ns, fittedvaluesset, 3;
        datacolour, fittedcolour, markersize, hidex=true, xticklabelrotation, xticks, xtitle
    )
    plotrenewalequationsamples_causaleffect!(
        gl, cases, cases_counterfactual, Ns, fittedvaluesset, 4;
        simcolour, fittedcolour, markersize, rhoclip, 
        hidex=false, xticklabelrotation, xticks, xtitle
    )

    colgap!(gl, 1, 5)
    if isnothing(columntitles)
        rowgap!(gl, 4, 5)
    else
        for r ∈ [ 1, 5 ] rowgap!(gl, r, 5) end 
    end
end

function plotrenewalequationsamples_w!(
    gl::GridLayout, 
    cases::Matrix, w, fittedvaluesset, fittedws, row;
    datacolour=:black, fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3,
    hidex=true,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle=L"$w_{jt}$",
)
    nlocations = size(cases, 2)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]

    for i ∈ 1:nlocations
        ws = _modelquantiles(fittedws, i)     
        band!(axs[i], xs, ws[:, 1], ws[:, 3]; color=fittedbandcolour)
        lines!(axs[i], xs, ws[:, 2]; color=fittedcolour, linewidth=1,)
        scatter!(
            axs[i], [ RenewalDiffInDiff._skip(x) ? missing : x for x ∈ w[:, i] ];
            color=datacolour, markersize
        )
    end
    linkyaxes!(axs...)

    for col ∈ 1:nlocations
        if col == 1 
            formataxis!(axs[col]; hidex)
        else
            formataxis!(axs[col]; hidex, hidey=true,)
        end
    end

    if !isnothing(ytitle)
        Label(
            gl[row, 0], ytitle; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    end
    if !hidex && !isnothing(xtitle)
        Label(gl[(row + 1), 1:nlocations], xtitle; fontsize=11.84, tellwidth=false)
    end

    return axs
end

function plotrenewalequationsamples_r0!(
    gl::GridLayout, 
    cases::Matrix, fittedvaluesset, row;
    betafunctions=nothing, 
    simcolour=COLOURVECTOR[3], 
    fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    infectiousduration=nothing, markersize=3, rhoclip=Inf,
    hidex=true,
    locationinds=RenewalDiffInDiff.automatic,
    plotcounterfactuals=false, 
    counterfactualcolour=( COLOURVECTOR[2], 0.75 ), 
    counterfactualbandcolour=( counterfactualcolour, 0.5 ),
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle=L"\mathcal{R}_0",
    yticks=Makie.automatic, 
)
    if locationinds isa RenewalDiffInDiff.Automatic
        nlocations = size(cases, 2)
        _locinds = 1:nlocations
    else
        nlocations = length(locationinds)
        _locinds = locationinds
    end

    duration = size(cases, 1)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks, yticks) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_locinds)
        rhoqs = _modelquantiles(fittedvaluesset, :rho_matrix_vec, ℓ)
        inds = findall(x -> x <= rhoclip, rhoqs[:, 3])
        band!(axs[i], xs[inds], rhoqs[:, 1][inds], rhoqs[:, 3][inds]; color=fittedbandcolour)
        lines!(axs[i], xs[inds], rhoqs[:, 2][inds]; color=fittedcolour, linewidth=1,)
        if plotcounterfactuals 
            crhoqs = _modelquantiles(fittedvaluesset, :rho_matrix_vec_counterfactual, ℓ)
            cinds = findall(x -> x <= rhoclip, crhoqs[:, 2])
            band!(
                axs[i], xs[cinds], crhoqs[:, 1][cinds], crhoqs[:, 3][cinds]; 
                color=counterfactualbandcolour
            )  
            lines!(
                axs[i], xs[cinds], crhoqs[:, 2][cinds]; 
                color=counterfactualcolour, linewidth=1,
            )    
        end

        isnothing(betafunctions) && continue
        scatter!(
            axs[i], [ betafunctions[ℓ](t) * infectiousduration for t ∈ 1:duration ]; 
            color=simcolour, markersize
        )
    end
    linkyaxes!(axs...)

    for col ∈ 1:nlocations
        if col == 1 
            formataxis!(axs[col]; hidex,)
        else
            formataxis!(axs[col]; hidex, hidey=true,)
        end
    end

    if !isnothing(ytitle)
        Label(
            gl[row, 0], ytitle; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    end

    if !hidex && !isnothing(xtitle)
        Label(gl[(row + 1), 1:nlocations], xtitle; fontsize=11.84, tellwidth=false)
    end

    return axs
end

function plotrenewalequationsamples_cases!(
    gl::GridLayout, 
    cases::Matrix, Ns::Vector, fittedvaluesset, row;
    fittedparameter=:y_matrix_poisson_vec,
    datacolour=:black, fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3, 
    hidex=true, 
    locationinds=RenewalDiffInDiff.automatic,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle="Diagnoses",
)
    if locationinds isa RenewalDiffInDiff.Automatic
        nlocations = size(cases, 2)
        _locinds = 1:nlocations
    else
        nlocations = length(locationinds)
        _locinds = locationinds
    end

    duration = size(cases, 1)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_locinds)
        yqs = _modelquantiles(fittedvaluesset, fittedparameter, ℓ)

        band!(
            axs[i], xs, 100_000 .* yqs[:, 1] ./ Ns[ℓ], 100_000 .* yqs[:, 3] ./ Ns[ℓ]; 
            color=fittedbandcolour
        )
        lines!(
            axs[i], xs, 100_000 .* yqs[:, 2] ./ Ns[ℓ]; 
            color=fittedcolour, linewidth=1,
        )
        scatter!(
            axs[i], 100_000 .* cases[:, ℓ] ./ Ns[ℓ]; 
            color=datacolour, markersize
        )
    end
    linkyaxes!(axs...)

    for col ∈ 1:nlocations
        if col == 1 
            formataxis!(axs[col]; hidex,)
        else
            formataxis!(axs[col]; hidex, hidey=true,)
        end
    end

    if !isnothing(ytitle)
        Label(
            gl[row, 0], ytitle; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    end

    if !hidex && !isnothing(xtitle)
        Label(gl[(row + 1), 1:nlocations], xtitle; fontsize=11.84, tellwidth=false)
    end

    return axs
end

function plotrenewalequationsamples_causaleffect!(
    gl::GridLayout, 
    cases::Matrix, cases_counterfactual, Ns::Vector, fittedvaluesset, row;
    fittedparameter=:y_matrix_poisson_vec,
    counterfactualfittedparameter=:y_matrix_poisson_vec_counterfactual,
    cumulativedifference=false,
    simcolour=COLOURVECTOR[3], 
    fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3, rhoclip=Inf,
    locationinds=RenewalDiffInDiff.automatic,
    hidex=false,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle="Difference",
    ytickformat=Makie.automatic,
)
    if locationinds isa RenewalDiffInDiff.Automatic
        nlocations = size(cases, 2)
        _locinds = 1:nlocations
    else
        nlocations = length(locationinds)
        _locinds = locationinds
    end

    duration = size(cases, 1)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks, ytickformat) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_locinds)
        yqs_cf = _modelquantiles(
            fittedvaluesset, fittedparameter, counterfactualfittedparameter, ℓ
        )
        if cumulativedifference 
            band!(
                axs[i], 
                xs, 
                100_000 .* cumsum(yqs_cf[:, 1]) ./ Ns[ℓ], 
                100_000 .* cumsum(yqs_cf[:, 3]) ./ Ns[ℓ]; 
                color=fittedbandcolour
            )
            lines!(
                axs[i], xs, 100_000 .* cumsum(yqs_cf[:, 2]) ./ Ns[ℓ]; 
                color=fittedcolour, linewidth=1,
            )
        else 
            band!(
                axs[i], 
                xs, 
                100_000 .* yqs_cf[:, 1] ./ Ns[ℓ], 
                100_000 .* yqs_cf[:, 3] ./ Ns[ℓ]; 
                color=fittedbandcolour
            )
            lines!(
                axs[i], xs, 100_000 .* yqs_cf[:, 2] ./ Ns[ℓ]; 
                color=fittedcolour, linewidth=1,
            )
        end
        
        isnothing(cases_counterfactual) && continue
        if cumulativedifference 
            scatter!(
                axs[i], 
                100_000 .* cumsum(cases[:, ℓ] .- cases_counterfactual[:, ℓ]) ./ Ns[ℓ]; 
                color=simcolour, markersize
            )
        else 
            scatter!(
                axs[i], 
                100_000 .* (cases[:, ℓ] .- cases_counterfactual[:, ℓ]) ./ Ns[ℓ]; 
                color=simcolour, markersize
            )
        end
    end
    linkyaxes!(axs...)

    for col ∈ 1:nlocations
        if col == 1 
            formataxis!(axs[col]; hidex,)
        else
            formataxis!(axs[col]; hidex, hidey=true,)
        end
    end

    if !isnothing(ytitle)
        Label(
            gl[row, 0], ytitle; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    end
    if !hidex && !isnothing(xtitle)
        Label(gl[(row + 1), 1:nlocations], xtitle; fontsize=11.84, tellwidth=false)
    end

    return axs
end

function _modelquantiles(vec, col; quantiles=[ 0.05, 0.5, 0.95 ])
    output = Matrix{Float64}(undef, length(vec[1][:, col]), length(quantiles))
    for i ∈ axes(output, 1)
        output[i, :] = quantile(
            [ vec[x][i, col] for x ∈ eachindex(vec) ],
            quantiles
        )
    end
    return output
end

function _modelquantiles(dataset, var, col; quantiles=[ 0.05, 0.5, 0.95 ])
    output = Matrix{Float64}(
        undef, length(getproperty(dataset, var)[1][:, col]), length(quantiles)
    )
    for i ∈ axes(output, 1)
        output[i, :] = quantile(
            [ 
                getproperty(dataset, var)[x][i, col] 
                for x ∈ eachindex(getproperty(dataset, var)) 
            ],
            quantiles
        )
    end
    return output
end

function _modelquantiles(dataset, var, var_counterfactual, col; quantiles=[ 0.05, 0.5, 0.95 ])
    output = Matrix{Float64}(
        undef, length(getproperty(dataset, var)[1][:, col]), length(quantiles)
    )
    for i ∈ axes(output, 1)
        output[i, :] = quantile(
            [ 
                getproperty(dataset, var)[x][i, col] - 
                    getproperty(dataset, var_counterfactual)[x][i, col] 
                for x ∈ eachindex(getproperty(dataset, var)) 
            ],
            quantiles
        )
    end
    return output
end

function fitws(cases::Matrix, Ns::Vector, fittedvaluesset)
    @unpack psi_vec, rho_matrix_vec = fittedvaluesset 
    ws1 = Matrix{Float64}(undef, size(rho_matrix_vec[1]))
    for t ∈ axes(ws1, 1), g ∈ axes(ws1, 2)
        if t == 1 
            ws1[t, g] = log(rho_matrix_vec[1][t, g]) + log(1)
        else
            ws1[t, g] = log(rho_matrix_vec[1][t, g]) + 
                log(1 - sum(cases[1:(t - 1), g]) / (psi_vec[1] * Ns[g]))
        end
    end
    w_vec = Vector{typeof(ws1)}(undef, length(psi_vec))
    for j ∈ eachindex(psi_vec)
        if j == 1 
            w_vec[j] = deepcopy(ws1)
        else 
            for t ∈ axes(ws1, 1), g ∈ axes(ws1, 2)
                if t == 1 
                    ws1[t, g] = log(rho_matrix_vec[j][t, g]) + log(1)
                else
                    ws1[t, g] =+(
                        log(rho_matrix_vec[j][t, g]),
                        log(1 - sum(cases[1:(t - 1), g]) / (psi_vec[j] * Ns[g]))  # this should already be constrained to be positive
                    )  
                        
                end
            end
            w_vec[j] = deepcopy(ws1)
        end
    end
    return w_vec
end
