
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

function plotchains(
    data::DataFrame; 
    size=( 400, 5000 ), ylabels=RenewalDiffInDiff.automatic, kwargs...
)
    @unpack colnames, plotnames_ind = _processplotchains(data; kwargs...)
        
    fig = Figure(; size)
    ax = [ Axis(fig[i, 1]) for i ∈ eachindex(plotnames_ind) ]
    n_ax = length(ax)
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j])
            if j == 1 
                if ylabels == RenewalDiffInDiff.automatic
                    Label(fig.layout[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
                else 
                    Label(fig.layout[i, 0], ylabels[i]; tellheight=false)
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

    Label(fig[n_ax+1, 1], "sample"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, n_ax, 5)
    
    return fig
end

function _processplotchains(data; logdensity="log_density")
    # "log_density" is the label given by `Pigeons` output. Turing labels it "lp".
    colnames = names(data)
    lp_ind = findall(x -> x == logdensity, colnames)
    _plotnames_ind = findall(x -> x ∉ _NOPLOTNAMES, colnames)
    plotnames_ind = [ lp_ind; _plotnames_ind ]
    return @ntuple colnames plotnames_ind
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
    infectiousduration=1, markersize=3, rhoclip=Inf,
    columntitles=nothing, columntitlefontsize=11.84, 
    xticklabelrotation=0.0, xticks=Makie.automatic, xtitle="Time"
)
    duration = size(cases, 1)
    nlocations = size(cases, 2)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])

    axs1 = [ Axis(gl[1, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]
    axs2 = [ Axis(gl[2, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]
    axs3 = [ Axis(gl[3, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]
    axs4 = [ Axis(gl[4, i]; xticklabelrotation, xticks) for i ∈ 1:nlocations ]

    for i ∈ 1:nlocations
        ws = _modelquantiles(fittedws, i)
        rhoqs = _modelquantiles(fittedvaluesset, :rho_matrix_vec, i)
        yqs = _modelquantiles(fittedvaluesset, :y_matrix_poisson_vec, i)
        yqs_cf = _modelquantiles(
            fittedvaluesset, :y_matrix_poisson_vec, :y_matrix_poisson_vec_counterfactual, i
        )

        if !isnothing(columntitles)
            Label(
            gl[0, i], columntitles[i]; 
            fontsize=columntitlefontsize, halign=:left, tellwidth=false
        )
        end
        
        band!(axs1[i], xs, ws[:, 1], ws[:, 2]; color=fittedcolour)
        scatter!(
            axs1[i], [ RenewalDiffInDiff._skip(x) ? missing : x for x ∈ w[:, i] ];
            color=datacolour, markersize
        )

        inds = findall(x -> x <= rhoclip, rhoqs[:, 2])
        band!(axs2[i], xs[inds], rhoqs[:, 1][inds], rhoqs[:, 2][inds]; color=fittedcolour)

        band!(
            axs3[i], xs, 100_000 .* yqs[:, 1] ./ Ns[i], 100_000 .* yqs[:, 2] ./ Ns[i]; 
            color=fittedcolour
        )
        scatter!(
            axs3[i], 100_000 .* cases[:, i] ./ Ns[i]; 
            color=datacolour, markersize
        )

        band!(
            axs4[i], xs, 100_000 .* yqs_cf[:, 1] ./ Ns[i], 100_000 .* yqs_cf[:, 2] ./ Ns[i]; 
            color=fittedcolour
        )
        
        isnothing(betafunctions) && continue
        scatter!(
            axs2[i], [ betafunctions[i](t) * infectiousduration for t ∈ 1:duration ]; 
            color=simcolour, markersize
        )
    
        isnothing(cases_counterfactual) && continue
        scatter!(
            axs4[i], 100_000 .* (cases[:, i] .- cases_counterfactual[:, i]) ./ Ns[i]; 
            color=simcolour, markersize
        )
    end
    linkyaxes!(axs1...)
    linkyaxes!(axs2...)
    linkyaxes!(axs3...)
    linkyaxes!(axs4...)

    for col ∈ 1:nlocations
        if col == 1 
            formataxis!(axs1[col]; hidex=true,)
            formataxis!(axs2[col]; hidex=true,)
            formataxis!(axs3[col]; hidex=true,)
            formataxis!(axs4[col])
        else
            formataxis!(axs1[col]; hidex=true, hidey=true,)
            formataxis!(axs2[col]; hidex=true, hidey=true,)
            formataxis!(axs3[col]; hidex=true, hidey=true,)
            formataxis!(axs4[col]; hidey=true,)
        end
    end

    Label(gl[1, 0], L"$w_{gt}$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(
        gl[2, 0], L"\mathcal{R}_0"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(gl[3, 0], "Diagnoses"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(
        gl[4, 0], "Difference"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(gl[5, 1:nlocations], xtitle; fontsize=11.84, tellwidth=false)

    colgap!(gl, 1, 5)
    if isnothing(columntitles)
        rowgap!(gl, 4, 5)
    else
        for r ∈ [ 1, 5 ] rowgap!(gl, r, 5) end 
    end
end

function _modelquantiles(vec, col; quantiles=[ 0.05, 0.95 ])
    output = Matrix{Float64}(undef, length(vec[1][:, col]), length(quantiles))
    for i ∈ axes(output, 1)
        output[i, :] = quantile(
            [ vec[x][i, col] for x ∈ eachindex(vec) ],
            quantiles
        )
    end
    return output
end

function _modelquantiles(dataset, var, col; quantiles=[ 0.05, 0.95 ])
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

function _modelquantiles(dataset, var, var_counterfactual, col; quantiles=[ 0.05, 0.95 ])
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
