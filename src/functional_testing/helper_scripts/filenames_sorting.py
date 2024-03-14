import editdistance
import os

def sort_filenames(filenames):
    """
    Enables users to just dump all files and this sorts these. So if there are forward and backward files for one sample, than this function will put both in a list.
    The function works based on levenshtein distance and if the distance is 1 between two basenames and the the distance is 0 if you cut off the ending of a filename (which is for fastq files 1 or 2) than these files will be matched
    :return: dict with matched forward and backward files and single strings for single files
    """ 
    filepath_basename = {}
    collected_basenames = set()
    for file in filenames:
        basename = os.path.basename(file).split('.')[0]

        if filepath_basename.get(basename) == None and basename not in collected_basenames:  
            filepath_basename[basename] = [file]
            collected_basenames.add(basename)
        
        for file2 in filenames:
            basename_keys = os.path.basename(file2).split('.')[0]
            ls_dist = editdistance.distance(basename, basename_keys)
            if ls_dist == 1:
                if (basename.endswith("_1") and basename_keys.endswith("_2")) or (basename.endswith("_2") and basename_keys.endswith("_1")):
                    cur_filepath = filepath_basename.get(basename)
                    if cur_filepath != None:
                        cur_filepath = cur_filepath[0]
                        matched_paths = [file, file2]
                        collected_basenames
                        if not basename_keys in collected_basenames:
                            filepath_basename[basename] = matched_paths
                            collected_basenames.add(basename)
                            collected_basenames.add(basename_keys)
                else:
                    pass                          
            else:
                pass

    return filepath_basename