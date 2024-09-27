


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# colour scheme
const COLOURVECTOR = [ 
    :blue, :seagreen4, :plum, :brown2, :darkgoldenrod1, :dodgerblue3, :skyblue2, :lightgray 
]

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
    @unpack colnames, plotnames_ind = _processplotchains(data; kwargs...)
        
    fig = Figure(; size)
    ax = [ Axis(fig[i, 1]) for i ∈ eachindex(plotnames_ind) ]
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j])
            Label(fig.layout[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
        end
    end
    
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

function plotrenewalequationsamples(dataset::Dict, w, fittedvaluesset; kwargs...)
    @unpack cases, cases_counterfactual, Ns = dataset 
    fittedws = fitws(cases, Ns, fittedvaluesset)
    return plotrenewalequationsamples(
        cases, cases_counterfactual, w, Ns, fittedvaluesset, fittedws; 
        kwargs...
    )
end

function plotrenewalequationsamples(cases::Matrix, w, Ns::Vector, fittedvaluesset; kwargs...)
    fittedws = fitws(cases, Ns, fittedvaluesset)
    return plotrenewalequationsamples(cases, w, Ns, fittedvaluesset, fittedws; kwargs...)
end

function plotrenewalequationsamples(
    cases::Matrix, w, Ns::Vector, fittedvaluesset, fittedws; 
    kwargs...
)
    return plotrenewalequationsamples(
        cases, nothing, w, Ns, fittedvaluesset, fittedws; 
        kwargs...
    )
end

function plotrenewalequationsamples(
    cases::Matrix, cases_counterfactual, w, Ns::Vector, fittedvaluesset, fittedws;
    betafunctions=nothing, betafunctions_counterfactual=nothing,
    datacolour=COLOURVECTOR[1], simcolour=COLOURVECTOR[2], fittedcolour=( :gray, 0.75 ), 
    infectiousduration=1, markersize=3, plotsize=( 800, 800 ),
)
    duration = size(cases, 1)
    nlocations = size(cases, 2)
    xs = eachindex(fittedvaluesset.rho_matrix_vec[1][:, 1])
    fig = Figure(; size=plotsize)
    axs1 = [ Axis(fig[1, i]) for i ∈ 1:nlocations ]
    axs2 = [ Axis(fig[2, i]) for i ∈ 1:nlocations ]
    axs3 = [ Axis(fig[3, i]) for i ∈ 1:nlocations ]
    axs4 = [ Axis(fig[4, i]) for i ∈ 1:nlocations ]
    #axs5 = [ Axis(fig[5, i]) for i ∈ 1:nlocations ]

    for i ∈ 1:nlocations
        ws = _modelquantiles(fittedws, i)
        rhoqs = _modelquantiles(fittedvaluesset, :rho_matrix_vec, i)
        yqs = _modelquantiles(fittedvaluesset, :y_matrix_poisson_vec, i)
        yqs_cf = _modelquantiles(
            fittedvaluesset, :y_matrix_poisson_vec, :y_matrix_poisson_vec_counterfactual, i
        )
        #logyqs_cf = _logmodelquantiles(fittedvaluesset, :y_matrix_poisson_vec, :y_matrix_poisson_vec_counterfactual, i)
        
        band!(axs1[i], xs, ws[:, 1], ws[:, 2]; color=fittedcolour)
        #for j ∈ eachindex(fittedvaluesset.rho_matrix_vec)
        #    lines!(axs1[i], log.(fittedvaluesset.rho_matrix_vec[j][:, i]); color=fittedcolour)
        #end
        scatter!(
            axs1[i], [ RenewalDiffInDiff._skip(x) ? missing : x for x ∈ w[:, i] ];
            color=datacolour, markersize
        )
    #end
    #for i ∈ 1:nlocations
        #for j ∈ eachindex(fittedvaluesset.rho_matrix_vec)
        #    lines!(axs2[i], fittedvaluesset.rho_matrix_vec[j][:, i]; color=fittedcolour)
        #end
        band!(axs2[i], xs, (rhoqs[:, 1]), (rhoqs[:, 2]); color=fittedcolour)

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
        
        #band!(axs5[i], xs, logyqs_cf[:, 1], logyqs_cf[:, 2]; color=fittedcolour)

        isnothing(betafunctions) && continue
        scatter!(
            axs2[i], [ betafunctions[i](t) * infectiousduration for t ∈ 1:duration ]; 
            color=simcolour, markersize
        )
    #end
    #for i ∈ 1:nlocations
        #for j ∈ eachindex(fittedvaluesset.y_matrix_poisson_vec) 
        #    lines!(
        #        axs3[i], 100_000 .* fittedvaluesset.y_matrix_poisson_vec[j][:, i] ./ Ns[i]; 
                #axs3[i], 100_000 .* fittedvaluesset.y_matrix_det_vec[j][:, i] ./ Ns[i]; 
        #        color=fittedcolour
        #    )
        #end

    #end 
    #=
    axs3 = [ Axis(fig[3, i]) for i ∈ 1:nlocations ]
    for i ∈ 1:nlocations
        for j ∈ eachindex(fittedvaluesset.rho_matrix_vec) 
            lines!(
                axs3[i], 
                fittedvaluesset.rho_matrix_vec[j][:, i] ./ 
                    fittedvaluesset.rho_matrix_vec[j][:, i]; 
                color=fittedcolour
            )
        end
        isnothing(betafunctions) && continue
        scatter!(
            axs3[i], 
            [ 
                betafunctions[i](t) ./ betafunctions_counterfactual[i](t) 
                for t ∈ 1:duration 
            ]; 
            color=simcolour, markersize
        )
    end =#
    #for i ∈ 1:nlocations 
    #    for j ∈ eachindex(fittedvaluesset.rho_matrix_vec) 
    #        lines!(
    #            axs4[i], 
    #            100_000 .* 
    #                (fittedvaluesset.y_matrix_poisson_vec[j][:, i] .- 
    #                fittedvaluesset.y_matrix_poisson_vec_counterfactual[j][:, i]) ./     
    #                #(fittedvaluesset.y_matrix_det_vec[j][:, i] .- 
    #                #fittedvaluesset.y_matrix_det_vec_counterfactual[j][:, i]) ./ 
    #                Ns[i]; 
    #            color=fittedcolour
    #        )
    #    end
        isnothing(cases_counterfactual) && continue
        scatter!(
            axs4[i], 
            100_000 .* 
                (cases[:, i] .- cases_counterfactual[:, i]) ./ 
                Ns[i]; 
            color=simcolour, markersize
        )
    end
    linkyaxes!(axs1...)
    linkyaxes!(axs2...)
    linkyaxes!(axs3...)
    linkyaxes!(axs4...)

    return fig
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
                    ws1[t, g] = log(rho_matrix_vec[j][t, g]) + 
                        log(1 - sum(cases[1:(t - 1), g]) / (psi_vec[j] * Ns[g]))
                end
            end
            w_vec[j] = deepcopy(ws1)
        end
    end
    return w_vec
end
