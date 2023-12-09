import editdistance
import os

def sort_filenames(filenames):
    """
    Enables users to just dump all files and this sorts these. So if there are forward and backward files for one sample, than this function will put both in a list.
    The function works based on levenshtein distance and if the distance is 1 between two basenames and the the distance is 0 if you cut off the ending of a filename (which is for fastq files 1 or 2) than these files will be matched
    :return: list with matched forward and backward files and single strings for single files
    """ 
    filepath_basename = {}
    for file in filenames:
        basename = os.path.basename(file).split('.')[0]
        keys_dic = list(filepath_basename.keys())
        if len(keys_dic) != 0:
            for basename_keys in keys_dic:
                ls_dist = editdistance.distance(basename, basename_keys)
                if ls_dist == 1:
                    if basename.endswith("_1") and basename_keys.endswith("_2"):
                        cur_filepath = filepath_basename.get(basename_keys)
                        cur_filepath = cur_filepath[0]
                        matched_paths = [cur_filepath, file]
                        filepath_basename[basename] = matched_paths
                    else:
                        pass                          
                else:
                    pass
            if filepath_basename.get(basename) == None:       
                filepath_basename[basename] = [file]
            
        else:
            filepath_basename[basename] = [file]
    return filepath_basename