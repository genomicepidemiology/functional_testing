import json

def check_output(output_json_base):
    json_file = output_json_base
    with open(json_file, 'r') as f:
        json_data_true = json.load(f)
    length_df =  len(json_data_true['seq_regions'])
    return length_df