""" 
    formataxis!(axis::Axis; <keyword arguments>)
    formataxis!(cb::Colorbar; horizontal=false)
    formataxis!(label::Label)
    formataxis!(legend::Legend; horizontal=true)
    formataxis!(axes::Array; <keyword arguments>)

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
    axis::Axis; 
    hidex=false, hidexticks=false, hidey=false, hideyticks=false, 
    setorigin=false, setpoint=nothing, 
    hidespines=nothing, trimspines=false, 
    ticksize=3.0,
    xgridvisible=false, ygridvisible=false,
)
    _formataxishidespines!(axis, hidespines)
    axis.spinewidth = 1
    axis.xtrimspine = trimspines
    axis.ytrimspine = trimspines
    axis.xgridstyle = (:dot, :dense)
    axis.xgridvisible = xgridvisible
    axis.xgridwidth = 1 
    axis.ygridstyle = (:dot, :dense)
    axis.ygridvisible = ygridvisible
    axis.ygridwidth = 1
    axis.xticksize = ticksize
    axis.yticksize = ticksize
    axis.xtickwidth = 1
    axis.ytickwidth = 1
    axis.xlabelsize = 12
    axis.ylabelsize = 12
    axis.xticklabelsize = 10
    axis.yticklabelsize = 10
    axis.titlealign = :left
    axis.titlesize = 12
    setorigin && setorigin!(axis) 
    setvalue!(axis, setpoint)
    hidex && hidexdecorations!(axis; ticks=hidexticks) 
    hidey && hideydecorations!(axis; ticks=hideyticks)  
    return axis
end 

function formataxis!(legend::Legend; horizontal=true, titleposition=RenewalDiD.automatic)
    legend.framevisible = false
    legend.titleposition = _legendtitleposition(titleposition, horizontal)
    legend.labelsize = 10
    legend.titlesize = 10
    legend.patchsize = (20, 20)
    if horizontal
        legend.orientation = :horizontal
    else 
        legend.margin = (10, 10, 10, 10)
    end 
    return legend
end 

function formataxis!(cb::Colorbar; horizontal=false)
    cb.ticklabelsize = 10
    cb.labelsize = 12
    if horizontal 
        cb.height = 10 
    else 
        cb.width = 10 
    end
    return cb
end 

formataxis!(label::Label) = label.fontsize = 12 

function formataxis!(axes::AbstractArray; kwargs...)
    for ax ∈ axes 
        formataxis!(ax; kwargs...) 
    end 
    return nothing
end 

_legendtitleposition(titleposition::Symbol, ::Any) = titleposition 
_legendtitleposition(::RenewalDiD.Automatic, horizontal) = horizontal ? :left : :top

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
setvalue!(axis::Axis, x, y) = scatter!(axis, [x], [y]; markersize=0)
setvalue!(::Any, ::Nothing) = nothing 
setvalue!(axis, xy::Tuple) = setvalue!(axis, xy...)
setorigin!(axis) = setvalue!(axis)

# hide spines
_formataxishidespines!(::Any, ::Nothing) = nothing
_formataxishidespines!(axis, hidespines::Symbol) = hidespines!(axis, hidespines)

function _formataxishidespines!(axis, hidespines::T) where T <: Union{<:AbstractArray, <:Tuple}
    for d ∈ hidespines 
        hidespines!(axis, d) 
    end
    return nothing
end 

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

function _labelplots!(
    labels::Vector{String}, layouts, rows::T, cols; 
    kwargs...
) where T <: Integer
    rowvector=(zeros(T, length(labels)) .+ rows) 
    return _labelplots!(labels, layouts, rowvector, cols; kwargs...)
end 

function _labelplots!(
    labels::Vector{String}, layouts, rows::Vector{<:Integer}, cols::T; 
    kwargs...
) where T <: Integer
    colvector = zeros(T, length(labels)) .+ cols 
    return _labelplots!(labels, layouts, rows, colvector; kwargs...)
end 

function _labelplots!(
    labels::Vector{String}, layouts::Vector, rows::Vector{<:Integer}, cols::Vector{<:Integer};
    kwargs...
)
    @assert length(labels) == length(layouts)
    for (row, col, label, layout) in zip(rows, cols, labels, layouts)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
    return nothing
end

function _labelplots!(
    labels::Vector{String}, layout, rows::Vector{<:Integer}, cols::Vector{<:Integer};
    kwargs...
)
    for (row, col, label) in zip(rows, cols, labels)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
    return nothing
end

function _labelplots!(
    label::String, layout, row::Integer, col::Integer;
    font="TeX Gyre Heros Bold", fontsize=14, halign=:left, padding=(0, 5, 5, 0)
)
    return Label(layout[row, col, TopLeft()], label; font, fontsize, halign, padding)
end 
