
using CairoMakie, DataFrames, Pigeons, PlotFormatting, Turing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const _NOPLOTNAMES = [  # outputs from MCMC that are not plotted 
    "iteration", 
    "chain", 
    "lp", 
    "n_steps", 
    "is_accept", 
    "acceptance_rate", 
    "log_density", 
    "hamiltonian_energy", 
    "hamiltonian_energy_error", 
    "max_hamiltonian_energy_error", 
    "tree_depth", 
    "numerical_error", 
    "step_size", 
    "nom_step_size"
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
    colnames=RenewalDiffInDiff.automatic, 
    plotnames_ind=RenewalDiffInDiff.automatic,
    ylabels=RenewalDiffInDiff.automatic, 
    ylabelrotation=π/2,
    kwargs...
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
            lines!(
                ax[i], getproperty(data, colnames[k])[inds]; 
                color=COLOURVECTOR[j], linewidth=1,
            )
            if j == 1 
                if ylabels == RenewalDiffInDiff.automatic
                    Label(
                        gl[i, 0], "$(colnames[k])"; 
                        rotation=ylabelrotation, tellheight=false,
                    )
                else 
                    Label(
                        gl[i, 0], ylabels[i]; 
                        rotation=ylabelrotation, tellheight=false,
                    )
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
    @unpack plotnames_ind = _processplotchains(data; logdensity)
    return @ntuple colnames=cn plotnames_ind
end

function _processplotchains(
    data, ::RenewalDiffInDiff.Automatic, pni; 
    logdensity="log_density"
)
    @unpack colnames = _processplotchains(data; logdensity)
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

function processplotchains(
    data; 
    colnames=RenewalDiffInDiff.automatic, 
    plotnames_ind=RenewalDiffInDiff.automatic, 
    logdensity="log_density"
)
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

function _plotlocinds(::RenewalDiffInDiff.Automatic, cases)
    nlocations = _plotnlocations(RenewalDiffInDiff.automatic, cases)
    return 1:nlocations
end 

_plotlocinds(locationinds, ::Any) = locationinds
_plotnlocations(::RenewalDiffInDiff.Automatic, cases) = size(cases, 2)
_plotnlocations(locationinds, ::Any) = length(locationinds)

function plotrenewalequationsamples_w!(
    gl::GridLayout, 
    cases::Matrix, w, fittedvaluesset, fittedws, row;
    locationinds=RenewalDiffInDiff.automatic,
    datacolour=:black, fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3,
    hidex=true,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle=L"$w_{jt}$",
    yticks=Makie.automatic, 
)
    nlocations = _plotnlocations(locationinds, cases)
    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks, yticks) for i ∈ 1:nlocations ]

    plotrenewalequationsamples_w!(
        axs, cases, w, fittedvaluesset, fittedws;
        locationinds, datacolour, fittedcolour, fittedbandcolour, markersize,
    )

    for (i, ℓ) ∈ enumerate(_plotlocinds(locationinds, cases))
        if i == 1 
            formataxis!(axs[i]; hidex)
        else
            formataxis!(axs[i]; hidex, hidey=true,)
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

function plotrenewalequationsamples_w!(
    axs::VecOrMat{<:Axis}, 
    cases::Matrix, w, fittedvaluesset, fittedws;
    locationinds=RenewalDiffInDiff.automatic,
    datacolour=:black, fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3,
    kwargs...
)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    for (i, ℓ) ∈ enumerate(_plotlocinds(locationinds, cases))
        ws = _modelquantiles(fittedws, ℓ)     
        band!(axs[i], xs, ws[:, 1], ws[:, 3]; color=fittedbandcolour)
        lines!(axs[i], xs, ws[:, 2]; color=fittedcolour, linewidth=1,)
        scatter!(
            axs[i], [ RenewalDiffInDiff._skip(x) ? missing : x for x ∈ w[:, ℓ] ];
            color=datacolour, markersize
        )
    end
    linkyaxes!(axs...)
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
    else
        nlocations = length(locationinds)
    end

    duration = size(cases, 1)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks, yticks) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_plotlocinds(locationinds, cases))
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
    counterfactualcases=nothing, simcolour=COLOURVECTOR[3], 
    fittedparameter=:y_matrix_poisson_vec,
    datacolour=:black, fittedcolour=COLOURVECTOR[1], fittedbandcolour=( fittedcolour, 0.5 ),
    markersize=3, 
    hidex=true, 
    locationinds=RenewalDiffInDiff.automatic,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle="Diagnoses",
)
    if locationinds isa RenewalDiffInDiff.Automatic
        nlocations = size(cases, 2)
    else
        nlocations = length(locationinds)
    end

    duration = size(cases, 1)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_plotlocinds(locationinds, cases))
        yqs = _modelquantiles(fittedvaluesset, fittedparameter, ℓ)

        band!(
            axs[i], xs, 100_000 .* yqs[:, 1] ./ Ns[ℓ], 100_000 .* yqs[:, 3] ./ Ns[ℓ]; 
            color=fittedbandcolour
        )
        lines!(
            axs[i], xs, 100_000 .* yqs[:, 2] ./ Ns[ℓ]; 
            color=fittedcolour, linewidth=1,
        )
        if !isnothing(counterfactualcases)
            scatter!(
                axs[i], 100_000 .* counterfactualcases[:, ℓ] ./ Ns[ℓ]; 
                color=simcolour, markersize
            )
        end
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
    locationinds=RenewalDiffInDiff.automatic,
    hidex=false,
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time", ytitle="Difference",
    ytickformat=Makie.automatic,
    kwargs...
)
    nlocations = _plotnlocations(locationinds, cases)
    duration = size(cases, 1)
    axs = [ Axis(gl[row, i]; xticklabelrotation, xticks, ytickformat) for i ∈ 1:nlocations ]

    for (i, ℓ) ∈ enumerate(_plotlocinds(locationinds, cases))
        plotrenewalequationsamples_causaleffect!(
            axs[i], 
            cases, cases_counterfactual, Ns, fittedvaluesset, ℓ;
            locationinds,
            hidex,
            xticklabelrotation, xticks, xtitle, ytitle,
            ytickformat,
            kwargs...
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
    ax::Axis, 
    cases::Matrix, cases_counterfactual, Ns::Vector, fittedvaluesset, index;
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
    yqs_cf = _modelquantiles(
        fittedvaluesset, fittedparameter, counterfactualfittedparameter, index;
        cumulativedifference
    )
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])
    
    band!(
        ax, 
        xs, 
        100_000 .* yqs_cf[:, 1] ./ Ns[index], 
        100_000 .* yqs_cf[:, 3] ./ Ns[index]; 
        color=fittedbandcolour
    )
    lines!(
        ax, xs, 100_000 .* yqs_cf[:, 2] ./ Ns[index]; 
        color=fittedcolour, linewidth=1,
    )
            
    _plotrenewalequationsamples_causaleffect_counterfactual!(
        ax, cases, cases_counterfactual, Ns, fittedvaluesset, index;
        cumulativedifference, simcolour, markersize,
    )
end

function _plotrenewalequationsamples_causaleffect_counterfactual!(
    ::Axis, ::Matrix, ::Nothing, ::Vector, ::Any, ::Any;
    kwargs...
)
    nothing
end

function _plotrenewalequationsamples_causaleffect_counterfactual!(
    ax::Axis, 
    cases::Matrix, cases_counterfactual, Ns::Vector, fittedvaluesset, index;
    cumulativedifference=false,
    simcolour=COLOURVECTOR[3], 
    markersize=3, 
)
    if cumulativedifference 
        scatter!(
            ax, 
            100_000 .* cumsum(cases[:, index] .- cases_counterfactual[:, index]) ./ Ns[index]; 
            color=simcolour, markersize
        )
    else 
        scatter!(
            ax, 
            100_000 .* (cases[:, index] .- cases_counterfactual[:, index]) ./ Ns[index]; 
            color=simcolour, markersize
        )
    end
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

function _modelquantiles(
    dataset, var, var_counterfactual, col; 
    quantiles=[ 0.05, 0.5, 0.95 ], cumulativedifference=false
)
    output = Matrix{Float64}(
        undef, length(getproperty(dataset, var)[1][:, col]), length(quantiles)
    )
    if cumulativedifference 
        for i ∈ axes(output, 1)
            output[i, :] = quantile(
                [ 
                    sum(@view getproperty(dataset, var)[x][1:i, col]) - 
                        sum(@view getproperty(dataset, var_counterfactual)[x][1:i, col]) 
                    for x ∈ eachindex(getproperty(dataset, var)) 
                ],
                quantiles
            )
        end
    else 
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

function plottwophi_att!(gl, df1, df2, leadlagtimes; kwargs...)
    plotphi_att!(gl, df1, df2, leadlagtimes; kwargs...)
end

function plottau_att!(gl::GridLayout, df::DataFrame, leadlagtimes; row=1, kwargs...)
    ax = Axis(gl[row, 1])
    _plotonephi_att!(ax, df, leadlagtimes; kwargs...)
    _labelplotphi_att!(gl, row, 1; kwargs...)
end

function plotphi_att!(gl::GridLayout, df1::DataFrame, df2::DataFrame, leadlagtimes; kwargs...)
    plotphi_att!(gl::GridLayout, [ df1, df2 ], leadlagtimes; kwargs...)
end

function plotphi_att!(gl::GridLayout, dfs::Vector{<:DataFrame}, leadlagtimes; row=1, kwargs...)
    axs = [ Axis(gl[row, i]) for i ∈ eachindex(dfs) ]
    for (i, df) ∈ enumerate(dfs)
        _plotonephi_att!(axs[i], df, leadlagtimes; hidey=(i != 1), kwargs...)
    end
    linkaxes!(axs...)
    _labelplotphi_att!(gl, row, 2; kwargs...)
end

function plotphi_att!(ax::Axis, df::DataFrame, leadlagtimes; kwargs...)
    _plotonephi_att!(ax, df, leadlagtimes; kwargs...)
end

function _plotonephi_att!(ax, df, leadlagtimes; hidey=false, kwargs...)
    tq = phiquantiles(df, leadlagtimes; kwargs...)
    scatter!(ax, leadlagtimes, tq.med; color=COLOURVECTOR[1], marker=:x, markersize=5,)
    rangebars!(ax, leadlagtimes, tq.lci, tq.uci; color=COLOURVECTOR[1], linewidth=1,)
    hlines!(
        ax, 0; 
        color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
    )
    formataxis!(
        ax; 
        hidespines=( :r, :t ), hidey, hideyticks=hidey, trimspines=true,
    )
    if hidey
        hidespines!(ax, :l)
    end
end

labelplotphi_att!(gl, row, width; kwargs...) = _labelplotphi_att!(gl, row, width; kwargs...)

function _labelplotphi_att!(gl, row, width; kwargs...)
    Label(
        gl[row, 0], L"$\varphi_{\mathrm{ATT}}$"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        gl[row+1, 1:width], "Time, relative to intervention"; 
        fontsize=11.84, tellwidth=false,
    )
    colgap!(gl, 1, 5)  
    rowgap!(gl, row, 5)
end