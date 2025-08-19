
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
    axis.xgridstyle = ( :dot, :dense )
    axis.xgridvisible = xgridvisible
    axis.xgridwidth = 1 
    axis.ygridstyle = ( :dot, :dense )
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
    if setorigin setorigin!(axis) end 
    setvalue!(axis, setpoint)
    if hidex hidexdecorations!(axis; ticks=hidexticks) end
    if hidey hideydecorations!(axis; ticks=hideyticks) end 
end 

function formataxis!(axis::Axis3; setorigin=false)
    axis.xspinewidth = 1 
    axis.yspinewidth = 1
    axis.zspinewidth = 1
    axis.xgridvisible = false
    axis.ygridvisible = false
    axis.zgridvisible = false
    axis.xtickwidth = 1
    axis.ytickwidth = 1
    axis.ztickwidth = 1
    axis.xlabelsize = 12
    axis.ylabelsize = 12
    axis.zlabelsize = 12
    axis.xticklabelsize = 10
    axis.yticklabelsize = 10
    axis.zticklabelsize = 10
    axis.titlealign = :left
    axis.titlesize = 12
    if setorigin setorigin!(axis) end
end 

function formataxis!(legend::Legend; horizontal=true, titleposition=automatic)
    legend.framevisible = false
    legend.labelsize = 10
    legend.titlesize = 10
    legend.patchsize = ( 20, 20 )
    if horizontal
        legend.orientation = :horizontal
        if titleposition == automatic
            legend.titleposition = :left
        else
            legend.titleposition = titleposition 
        end
    else 
        legend.margin = ( 10, 10, 10, 10 )
        if titleposition == automatic
            legend.titleposition = :top
        else
            legend.titleposition = titleposition 
        end
    end 
end 

function formataxis!(cb::Colorbar; horizontal=false)
    cb.ticklabelsize = 10
    cb.labelsize = 12
    if horizontal 
        cb.height = 10 
    else 
        cb.width = 10 
    end
end 

formataxis!(label::Label) = label.fontsize = 12 

function formataxis!(axes::AbstractArray; kwargs...)
    for ax ∈ axes 
        formataxis!(ax; kwargs...) 
    end 
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
setvalue!(axis::Axis, x, y) = scatter!(axis, [ x ], [ y ]; markersize=0)
setvalue!(axis::Axis3, z::Real=0) = setvalue!(axis, 0, 0, z)
setvalue!(axis::Axis3, x, y, z) = scatter!(axis, [ x ], [ y ], [ z ]; markersize=0)
setvalue!(::Any, ::Nothing) = nothing 
setvalue!(axis, xy::Tuple) = setvalue!(axis, xy...)
setorigin!(axis) = setvalue!(axis)

# Function to hide spines. Not exported.

_formataxishidespines!(axis, hidespines::Nothing) = nothing
_formataxishidespines!(axis, hidespines::Symbol) = hidespines!(axis, hidespines)

function _formataxishidespines!(axis, hidespines::T) where T <: Union{<:AbstractArray, <:Tuple}
    for d ∈ hidespines hidespines!(axis, d) end
end 
