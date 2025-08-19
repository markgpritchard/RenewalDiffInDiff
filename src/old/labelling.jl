
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
