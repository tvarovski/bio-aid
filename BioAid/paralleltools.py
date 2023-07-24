# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import os
import argparse
import logging
import subprocess
import multiprocessing
from math import ceil
import pandas as pd

def splitDF(df: pd.DataFrame, slices_num: int, savecsv_path_root=None) -> list:
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

def joinSlices(slices_num: int, slice_path_root: str, savecsv_path=None) -> pd.DataFrame:
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

def cleanUpSlices(slices_num: int, slice_path_root: str) -> None:
    """
    Deletes the CSV files for each slice of a DataFrame that was split using the `split_df` function.

    Args:
        slices_num (int): The number of slices that were created.
        slice_path_root (str): The root path to the CSV files for each slice.

    Returns:
        None
    """
    logging.info('Combining results from each slice...')
    for slice in range(slices_num):
        os.remove(f"{slice_path_root}_{slice}.csv")
    logging.info(f'Clean up Finished! Deleted slice files.')

def runScriptSubprocess(args: tuple[str, str]) -> None:
    """
    Runs a Python script in a subprocess using the provided CSV file as input.

    Args:
        args (tuple): A tuple containing the path to the CSV file and the path to the Python script.

    Returns:
        None
    """
    csv_slice, script_path = args
    subprocess.call(['python', script_path, csv_slice])

def runMainPool(script_path: str, results_df_path: str, num_processes: int = 8, joint_name: str = 'all') -> None:
    """
    Runs a Python script in parallel on a pandas DataFrame using the provided CSV file as input.

    Args:
        script_path (str): The path to the Python script to be run.
        results_df_path (str): The path to the CSV file containing the input DataFrame.
        num_processes (int): Optional. The number of processes to use for parallelization. Default is 8.
        joint_name (str): Optional. The name of the output CSV file containing the combined results. Default is 'all'.

    Returns:
        None
    """
    with open(results_df_path, 'r') as f:
        df = pd.read_csv(f)
        results_df_name = os.path.basename(results_df_path)[:-4]

        # split the df into n=num_processes slices and save each slice to a csv file
        logging.info(f'Splitting {results_df_name} into {num_processes} slices...')
        splitDF(df, num_processes, results_df_name)

        # run the script in parallel on each slice
        logging.info(f'Running {script_path} in parallel on {num_processes} slices...')
        pool = multiprocessing.Pool(processes=num_processes)
        pool.map(runScriptSubprocess, [(f"{results_df_name}_{slice}.csv", script_path) for slice in range(num_processes)])
        pool.close()
        pool.join()
        logging.info(f'Finished {num_processes} processes...')

        # join the results from each slice into one DataFrame and save it to a CSV file
        joint_path = f'{results_df_name}_{joint_name}.csv'
        logging.info(f'Joining results from {num_processes} slices into {joint_path}...')
        joinSlices(num_processes, results_df_name, joint_path)

        # clean up the slice files
        logging.info(f'Cleaning up slice files for {results_df_name}...')
        cleanUpSlices(num_processes, results_df_name)

    logging.info(f'Finished processing {results_df_name}.')

def parseArguments() -> argparse.Namespace:
    """
    Parses command line arguments.

    Returns:
    Parsed arguments
    """
    program_description =  'This program runs a Python script in parallel on a pandas DataFrame using \
                            the provided CSV file as input. The CSV file must contain a header row. \
                            Caution: The CSV file will be split into n=num_processes slices and each slice \
                            will be saved to a CSV file. The slice files will be deleted after the script \
                            is finished. Make sure that your input file can be divided into sub-files \
                            without losing data...'
    
    parser = argparse.ArgumentParser(description=program_description)

    args_dict = {
        '-v': {'name': '--version', 'action': 'version', 'version': '%(prog)s 1.0'},
        '-t': {'name': '--target', 'type': str, 'default': 'myscript.py', 'help': 'Path to the Python script to be run'},
        '-r': {'name': '--results_df_path', 'type': str, 'default': 'results_df.csv', 'help': 'Path to the results dataframe to be processed'},
        '-p': {'name': '--num_processes', 'type': int, 'default': 8, 'help': 'Number of processes (cores) to use'}
    }

    for arg, properties in args_dict.items():
        parser.add_argument(arg, properties['name'], **{key: value for key, value in properties.items() if key != 'name'})

    return parser.parse_args()

def main() -> None:
    """
    This function is the entry point of the program. It sets up logging, and calls the function to run the main pool of processes.

    Returns:
    None
    """
    args = parseArguments()

    logging.info(f'Running {args.target} on {args.results_df_path} in parallel on {args.num_processes} processes...')
    runMainPool(args.target, args.results_df_path, num_processes=args.num_processes)
    logging.info(f'Finished running {args.target} on {args.results_df_path} in parallel on {args.num_processes} processes...')
    logging.info(f'Output saved to {args.results_df_path[:-4]}_all.csv')

if __name__ == "__main__":
    main()