
function addoffsetstointerventionarray(
    originalinterventions::InterventionMatrix{T}; 
    offset=-35:7:35
) where T 
    return InterventionArray(originalinterventions; offset)
end

function addoffsetstointerventionarray(
    originalinterventions::InterventionArray{T}; 
    offset=-35:7:35
) where T 
    _firstinterventions = InterventionMatrix{T}(
        originalinterventions.duration, originalinterventions.rawstarttimes[:, 1]
    )
    interventions = InterventionArray(_firstinterventions; offset)
    for k in axes(originalinterventions, 3)
        k == 1 && continue 
        interventions = cat(
            interventions, 
            InterventionMatrix{T}(
                originalinterventions.duration, originalinterventions.rawstarttimes[:, k]
            ); 
            dims=3
        )
    end
    return interventions
end
