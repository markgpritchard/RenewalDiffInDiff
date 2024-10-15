
using CairoMakie, DataFrames, Pigeons, Turing

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# colour scheme
const COLOURVECTOR = [
    RGBf(33 / 255, 145 / 255, 140 / 255),
    RGBf(68 / 255, 57 / 255, 131 / 255),
    RGBf(253 / 255, 231 / 255, 37 / 255),
    RGBf(53 / 255, 183 / 255, 121 / 255),
    RGBf(49 / 255, 104 / 255, 142 / 255),
    RGBf(68 / 255, 1 / 255, 84 / 255),
    RGBf(144 / 255, 215 / 255, 67 / 255),
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

function plotchains(data::DataFrame; size=( 400, 5000 ), ylabels=RenewalDiffInDiff.automatic, kwargs...)
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
            formataxis!(a; hidex=true, hidexticks=true, hidespines=( :b, :t, :r))
        end
    end

    Label(fig[n_ax+1, 1], "sample"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, n_ax, 5)

    #Label(fig[1, 0], "generation interval, f(τ)"; fontsize=11.84, rotation=π/2, tellheight=false)
    
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
    datacolour=COLOURVECTOR[1], simcolour=COLOURVECTOR[2], fittedcolour=( :gray, 0.75 ), 
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
            Label(gl[0, i], columntitles[i]; fontsize=columntitlefontsize, tellwidth=false)
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
            formataxis!(axs1[col]; hidex=true, hidexticks=true, hidespines=( :b, :r, :t ))
            formataxis!(axs2[col]; hidex=true, hidexticks=true, hidespines=( :b, :r, :t ))
            formataxis!(axs3[col]; hidex=true, hidexticks=true, hidespines=( :b, :r, :t ))
            formataxis!(axs4[col]; hidespines=( :r, :t ))
        else
            formataxis!(
                axs1[col]; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :b, :r, :t, :l )
            )
            formataxis!(
                axs2[col]; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :b, :r, :t, :l )
            )
            formataxis!(
                axs3[col]; 
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                hidespines=( :b, :r, :t, :l )
            )
            formataxis!(axs4[col]; hidey=true, hideyticks=true, hidespines=( :r, :t, :l ))
        end
    end

    Label(gl[1, 0], "w(gt)"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(
        gl[2, 0], "R₀"; 
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

""" 
    formataxis!(axis::Axis, width = 800; <keyword arguments>)
    formataxis!(axis::Axis3, width = 800; setorigin = false)
    formataxis!(cb::Colorbar, width = 800; horizontal = false)
    formataxis!(label::Label, width = 800)
    formataxis!(legend::Legend, width = 800; horizontal = true)
    formataxis!(axes::Array, args...; kwargs...)

Apply consistent formatting to components of figures.

## Arguments 
The first argument is the element to be formatted. 

`width` refers to the width of the whole figure. Text sizes and line widths are formatted 
    to be a constant proportion of this width.

## Keyword arguments
Most keyword arguments are only available when formatting an `Axis` 
* `hidespines = ( :r, :t )`: which sides of the plot are not to be displayed: accepts 
    a tuple of symbols, `:b` for bottom, `:l` for left, `:r` for right, `:t` for top
* `hidex = false`: whether to hide values on x-axis
* `hidexticks = false`: whether to hide tick marks on x-axis (can only be hidden if `hidex = true`)
* `hidey = false`: whether to hide values on y-axis
* `hideyticks = false`: whether to hide tick marks on y-axis (can only be hidden if `hidey = true`)
* `setorigin = false`: whether axes should be extended to a value of `0`
* `setpoint = nothing`: can take a number or tuple and calls `setvalue!`
* `trimspines = true`: whether to trim to axes at the minimum and maximum tick values
The additional keyword argument `horizontal` is available when formatting a `Colorbar`
    or `Legend`, and sets whether that item should be oriented horizontally.
""" 
function formataxis!(
    axis::Axis, width=800; 
    hidex=false, hidexticks=false, hidey=false, hideyticks=false, 
    setorigin=false, setpoint=nothing, #trimspines = true
)
    #formataxishidespines!(axis, hidespines)
    #axis.spinewidth = width / 800
    #axis.xtrimspine = trimspines; axis.ytrimspine = trimspines
    axis.xgridvisible = false; axis.ygridvisible = false
    axis.xtickwidth = width / 800; axis.ytickwidth = width / 800
    axis.xlabelsize = width / 67; axis.ylabelsize = width / 67
    axis.xticklabelsize = width / 80; axis.yticklabelsize = width / 80
    axis.titlealign = :left; axis.titlesize = width / 65
    if setorigin setorigin!(axis) end 
    setvalue!(axis, setpoint)
    if hidex 
        hidexdecorations!(axis; ticks = hidexticks) 
    else
        if hidexticks 
            @info "Function `formataxis!` cannot hide ticks on x axis unless `hidex` is true" 
        end 
    end 
    if hidey 
        hideydecorations!(axis; ticks = hideyticks) 
    else 
        if hideyticks 
            @info "Function `formataxis!` cannot hide ticks on y axis unless `hidey` is true" 
        end 
    end 
end 

function formataxis!(axis::Axis3, width = 800; setorigin = false)
    axis.xspinewidth = width / 800; axis.yspinewidth = width / 800; axis.zspinewidth = width / 800; 
    axis.xgridvisible = false; axis.ygridvisible = false; axis.zgridvisible = false
    axis.xtickwidth = width / 800; axis.ytickwidth = width / 800; axis.ztickwidth = width / 800
    axis.xlabelsize = width / 67; axis.ylabelsize = width / 67; axis.zlabelsize = width / 67
    axis.xticklabelsize = width / 80; axis.yticklabelsize = width / 80; axis.zticklabelsize = width / 80
    axis.titlealign = :left; axis.titlesize = width / 65
    if setorigin setorigin!(axis) end
end 

function formataxis!(legend::Legend, width = 800; horizontal = true)
    legend.framevisible = false
    legend.labelsize = width / 80; legend.titlesize = width / 80
    legend.patchsize = (width / 40, width / 40)
    if horizontal
        legend.orientation = :horizontal
        legend.titleposition = :left
    else 
        legend.margin = (10, 10, 10, 10)
    end 
end 

function formataxis!(cb::Colorbar, width = 800; horizontal = false)
    cb.ticklabelsize = width / 80; cb.labelsize = width / 67
    if horizontal 
        cb.height = width / 80 
    else 
        cb.width = width / 80 
    end
end 

function formataxis!(label::Label, width = 800)
    label.fontsize  = width / 67
end 

function formataxis!(axes::Array, width = 800; kwargs...)
    for ax ∈ axes formataxis!(ax, width; kwargs...) end 
end 

"""
    setvalue!(axis, <additional arguments>)

Extends axes to include the provided value. 

## Additional arguments 
Separate arguments can be given for the `x`, `y` and `z` (for `Axis3`) axes. They 
    may also be provided as a `Tuple`.  

For an `Axis3`, the `z` value may be provided alone, which assumes `x = y = 0`. 
    For an `Axis`, the `y` value may be provided alone, which assumes `x = 0`. If 
    no additional arguments are given, the origin, is added. If `x` is supplied as 
    `nothing`, no extension to the axes is performed.

"""
setvalue!(axis::Axis, y::Real=0) = setvalue!(axis, 0, y)
setvalue!(axis::Axis, x, y) = scatter!(axis, [x], [y], markersize = 0)
setvalue!(axis::Axis3, z::Real = 0) = setvalue!(axis, 0, 0, z)
setvalue!(axis::Axis3, x, y, z) = scatter!(axis, [x], [y], [z], markersize = 0)
setvalue!(::Any, ::Nothing) = nothing 
setvalue!(axis, xy::Tuple) = setvalue!(axis, xy...)
setorigin!(axis) = setvalue!(axis)

# Function to hide spines. Not exported.

formataxishidespines!(axis, hidespines::Nothing) = nothing
formataxishidespines!(axis, hidespines::Symbol) = hidespines!(axis, hidespines)

function formataxishidespines!(axis, hidespines)
    for d ∈ hidespines hidespines!(axis, d) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Label subplots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""" 
    labelplots!(labels, layouts; <keyword arguments>)

Applies labels to plots within a larger figure.

`labels` and `layouts` can each refer to an individual plot or can be a vector referring 
    to several plots.

# Keyword arguments 
* `cols = 0`: which column within the plots should have the label applied. Can be 
    provided as one integer for all plots or a vector of integers for each plot.
* `rows = 0`: which row within the plots should have the label applied. Can be provided 
    as one integer for all plots or a vector of integers for each plot.
* `font = "TeX Gyre Heros Bold"`: font for the labels
* `fontsize = 14`: label fontsize
* `halign = :left`: label alignment
* `padding = ( 0, 5, 5, 0 )`: label padding, provided as a `Tuple` in the order left, 
    right, bottom, top
"""
function labelplots!(labels, layouts; cols=0, rows=0, kwargs...)
    return _labelplots!(labels, layouts, rows, cols; kwargs...)
end 

function _labelplots!(labels::Vector{String}, layouts, rows::Int, cols; kwargs...)
    rowvector=(zeros(Int, length(labels)) .+ rows) 
    return _labelplots!(labels, layouts, rowvector, cols; kwargs...)
end 

function _labelplots!(
    labels::Vector{String}, layouts, rows::Vector{<:Int}, cols::Int; 
    kwargs...
)
    colvector = zeros(Int, length(labels)) .+ cols 
    return _labelplots!(labels, layouts, rows, colvector; kwargs...)
end 

function _labelplots!(
    labels::Vector{String}, layouts::Vector, rows::Vector{<:Int}, cols::Vector{<:Int};
    kwargs...
)
    @assert length(labels) == length(layouts)
    for (row, col, label, layout) ∈ zip(rows, cols, labels, layouts)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
end

function _labelplots!(
    labels::Vector{String}, layout, rows::Vector{<:Int}, cols::Vector{<:Int};
    kwargs...
)
    for (row, col, label) ∈ zip(rows, cols, labels)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
end

function _labelplots!(
    label::String, layout, row::Int, col::Int;
    font="TeX Gyre Heros Bold", fontsize=14, halign=:left, padding=( 0, 5, 5, 0 )
)
    Label(layout[row, col, TopLeft()], label; font, fontsize, halign, padding)
end 
