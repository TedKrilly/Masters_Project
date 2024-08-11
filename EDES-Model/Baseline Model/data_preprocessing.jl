using DataFrames

function preprocess_data(calibrated_params, complex_meal_data)
    common_ids = intersect(unique(calibrated_params.ID), unique(complex_meal_data.ID))
    
    filtered_params = filter(row -> row.ID in common_ids, calibrated_params)
    filtered_meal_data = filter(row -> row.ID in common_ids, complex_meal_data)
    
    valid_ids = []
    for id in common_ids
        id_data = filter(row -> row.ID == id, filtered_meal_data)
        if "gluc" in id_data.metab && "ins" in id_data.metab
            push!(valid_ids, id)
        end
    end
    
    final_params = filter(row -> row.ID in valid_ids, filtered_params)
    final_meal_data = filter(row -> row.ID in valid_ids, filtered_meal_data)
    
    sort!(final_meal_data, [:ID, :Condition, :time])
    
    return final_params, final_meal_data
end