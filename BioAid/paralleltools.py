# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## split_df
## join_slices
## clean_up_slices

import pandas as pd
import logging

def split_df(df: pd.DataFrame, slices_num: int, savecsv_path_root=None) -> list:
    """
    Splits a pandas DataFrame into n=slices_num slices and saves each slice to a CSV file if a path is provided.
    Returns a list of the DataFrame slices.

    Args:
        df (pandas.DataFrame): The DataFrame to be split.
        slices_num (int): The number of slices to create.
        savecsv_path_root (str): Optional. The root path to save the CSV files. If not provided, the slices are not saved.

    Returns:
        list: A list of the DataFrame slices.
    """

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

def join_slices(slices_num: int, slice_path_root: str, savecsv_path=None):
    """
    Combines the results from each slice into one pandas DataFrame and saves it to a CSV file if a path is provided.
    Returns the combined DataFrame.

    Args:
        slices_num (int): The number of slices to combine.
        slice_path_root (str): The root path to the CSV files for each slice.
        *savecsv_path (str): Optional. The path to save the combined DataFrame as a CSV file.

    Returns:
        pandas.DataFrame: The combined DataFrame.
    """

    df = pd.DataFrame()
    for slice in range(slices_num):
        df_slice = pd.read_csv(f"{slice_path_root}_{slice}.csv")
        df = pd.concat([df, df_slice], ignore_index=True)
    
    if savecsv_path:
        df.to_csv(savecsv_path[0], index=False)
        logging.info(f'Slices joined! Results saved to {savecsv_path[0]}')
    
    return df

def clean_up_slices(slices_num: int, slice_path_root: str) -> None:
    """
    Deletes the CSV files for each slice of a DataFrame that was split using the `split_df` function.

    Args:
        slices_num (int): The number of slices that were created.
        slice_path_root (str): The root path to the CSV files for each slice.

    Returns:
        None
    """
    import os
    for slice in range(slices_num):
        os.remove(f"{slice_path_root}_{slice}.csv")

    logging.info(f'Clean up Finished! Deleted slice files.')