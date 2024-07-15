from Exact.tolerances import *

def combine_exact_data(w_data,weight_data):
    combined_w_data = []
    combined_weight_data = []
    previous_w_value = None
    for index,value in enumerate(w_data):
        if previous_w_value == None:
            combined_w_data.append(value)
            combined_weight_data.append(weight_data[index])
            previous_w_value = value
        else:
            if abs(value - previous_w_value) <= TOLERANCE_combine:
                previous_index = len(combined_w_data) - 1
                combined_weight_data[previous_index] = combined_weight_data[previous_index] + weight_data[index]
                previous_w_value = value
            else:
                combined_w_data.append(value)
                combined_weight_data.append(weight_data[index])
                previous_w_value = value   
    return combined_w_data,combined_weight_data