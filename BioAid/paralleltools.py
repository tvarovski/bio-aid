# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## split_df
## join_slices
## clean_up_slices

import pandas as pd
import logging

def split_df(df, slices_num, *savecsv_path_root):
    #this function splits the df into n=slices_num slices and saves each slice to a csv file
    #if savecsv_path_root is given, the slices are saved to csv files with the given path root
    #and the slice number appended to the end of the path
    #the slices are returned as a list of df

    from math import ceil

    slice_size = ceil(len(df) / slices_num)
    slices_out = []
    for slice in range(slices_num):
        df_slice = df[slice*slice_size:(slice+1)*slice_size].copy()
        
        if savecsv_path_root:
            df_slice_name = f"{savecsv_path_root[0]}_{slice}.csv"
            df_slice.to_csv(df_slice_name, index=False)
            logging.info(f'Slice {slice} saved to {df_slice_name}')
        
        slices_out.append(df_slice)
    
    logging.info(f'Input df split into {slices_num} slices.')
    return slices_out

def join_slices(slices_num, slice_path_root, *savecsv_path):
    #this function combines the results from each slice into one df and saves it
    #if savecsv_path is given, the df is saved to the given path
    #the df is returned

    df = pd.DataFrame()
    for slice in range(slices_num):
        df_slice = pd.read_csv(f"{slice_path_root}_{slice}.csv")
        df = pd.concat([df, df_slice], ignore_index=True)
    
    if savecsv_path:
        df.to_csv(savecsv_path[0], index=False)
        logging.info(f'Slices joined! Results saved to {savecsv_path[0]}')
    
    return df

def clean_up_slices(slices_num, slice_path_root):
    import os
    #this function deletes the slice files

    for slice in range(slices_num):
        os.remove(f"{slice_path_root}_{slice}.csv")

    logging.info(f'Clean up Finished! Deleted slice files.')