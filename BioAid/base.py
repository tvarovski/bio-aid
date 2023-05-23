# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## extractSeqFromFastaToList
## validateSquence
## compl
## rev_compl
## unique
## createChrList
## dataFrameImport

import logging

def extractSeqFromFastaToList(fasta_file_path: str) -> list:
    '''Extracts sequence information from a FASTA file and returns it as a nested list.

    Args:
        fasta_file_path (str): The path to the FASTA file.

    Returns:
        list: A nested list of the form [[title_0, sequence_0], [title_1, sequence_1], ...].
            Each element of the list is a list containing the title and sequence of a sequence in the FASTA file.
    '''
    fasta_file = open(fasta_file_path, 'r')
    contents = fasta_file.readlines()
    fasta_file.close()

    fasta_list=[]

    for i in contents:
        if '>' in i:
            fasta_list.append([i.strip('>').strip(),''])
        else:
            fasta_list[-1][1] = fasta_list[-1][1]+i.strip()
    logging.info(f"Extraction of sequence information from {fasta_file_path} finished.")
    return fasta_list

def validateSquence(sequence: str) -> bool:
    '''Checks if a DNA sequence contains only valid nucleotides.

    Args:
        sequence (str): A DNA sequence to be validated.

    Returns:
        bool: True if the sequence contains only valid nucleotides, False otherwise.
            If the sequence contains non-canonical nucleotides, a log message is generated.
            If the sequence contains 'N's, a warning log message is generated.
            The log messages are written to the default logger.
    '''
    bases = "ATGCatgcN"
    for i in sequence:
        if i not in bases:
            logging.info(f"Sequence doesn't contain cannonical nucleotides: {i}")
            return(False)
        if i == "N":
            logging.warning(f"Warning, sequence contains 'N's")
    return(True)

def compl(base: str) -> str:
    """
    Returns the complementary base for a given DNA base.

    Args:
        base (str): A single character representing a DNA base. Must be one of 'A', 'T', 'G', 'C', or '-'.

    Returns:
        str: The complementary base for the given DNA base. Returns '-' if the input is '-'.

    Raises:
        ValueError: If the input is not one of 'A', 'T', 'G', 'C', or '-'.
    """
    if base == "A":
        return('T')
    elif base == "T":
        return('A')
    elif base == "G":
        return('C')
    elif base == "C":
        return('G')
    elif base == "-":
        return('-')
    else:
        logging.error(f"Invalid base: {base}")
        raise ValueError(f"Invalid base: {base}")

def rev_compl(seq: str) -> str:
    """
    Returns the reverse complement of a given DNA sequence.

    Args:
        seq (str): A string representing a DNA sequence. Must only contain characters 'A', 'T', 'G', 'C', or '-'.

    Returns:
        str: The reverse complement of the given DNA sequence.

    Raises:
        ValueError: If the input sequence contains characters other than 'A', 'T', 'G', 'C', or '-' (Through compl())
    """
    new_seq = ""
    for base in seq:
        new_base = compl(base)
        new_seq = new_base + new_seq
    return(new_seq)

def unique(list1: list) -> list:
    """
    Returns a new list containing only the unique elements of the input list.

    Args:
        list1 (list): A list of elements.

    Returns:
        list: A new list containing only the unique elements of the input list.

    Examples:
        >>> unique([1, 2, 3, 2, 1])
        [1, 2, 3]

        >>> unique(['a', 'b', 'c', 'b', 'a'])
        ['a', 'b', 'c']
    """

    list_set = set(list1) 
    unique_list = (list(list_set))
    return unique_list

def createChrList(chrNum: int) -> list:
    """
    Creates a list of chromosome names.

    Args:
        chrNum (int): The number of chromosomes to include in the list.

    Returns:
        list: A list of chromosome names, where each name is a string of the form "chrX", where X is the chromosome number.

    Examples:
        >>> createChrList(3)
        ['chr1', 'chr2', 'chr3']

        >>> createChrList(5)
        ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
    """
    chrList = []
    for i in range(chrNum):
        chrList.append(f"chr{i+1}")
    return(chrList)

def dataFrameImport(directory: str) -> tuple:
    """
    Imports all TSV files in a directory as Pandas dataframes.

    Args:
        directory (str): The path to the directory containing the TSV files.

    Returns:
        tuple: A tuple containing two elements:
            - A list of Pandas dataframes, where each dataframe corresponds to a TSV file in the directory.
            - A list of sample names, where each name is a string representing the name of the corresponding TSV file.
    """
    import os
    import pandas as pd
    frames_list_samples = []
    sample_names = []

    counter=1
    for filename in os.listdir(directory):
            if filename.endswith(".tsv"):
                    sample_names.append(filename[:10].strip())
                    path = os.path.join(directory, filename)
                    globals()[f"sample{counter}"] = pd.read_csv(path, sep='\t')
                    print(f'dataframe {sample_names[-1]} has {len(globals()[f"sample{counter}"])} total rows')
                    frames_list_samples.append(globals()[f"sample{counter}"])
                    counter+=1
    logging.info(f"found {len(frames_list_samples)} samples in {directory}")
    logging.info(len(sample_names), sample_names)

    return(frames_list_samples, sample_names)

def pullGenomicContext(list_of_positions: list, fasta_file_path: str, context_flank: int = 5) -> list:
    """
    Extracts genomic context of a specified length around each position in a list of positions.

    Args:
        list_of_positions (list): A list of positions, where each position is a list of two elements: the chromosome name and the position.
        fasta_file_path (str): The path to the FASTA file containing the reference genome.
        context_flank (int, optional): The length of the context to extract on either side of each position. Defaults to 5.

    Returns:
        list: A list of genomic contexts, where each context is a list of four elements:
            - The sequence of bases to the left of the position.
            - The query base at the position.
            - The sequence of bases to the right of the position.
            - The chromosome name.
    """
    fasta_list = extractSeqFromFastaToList(fasta_file_path)
    context_list = []

    for i in fasta_list:
        logging.debug(f"extracting context from {i[0]}. It has {len(i[1])} bases.")
        sequence = i[1]
        for j in list_of_positions:

            #check if it is the same chromosome
            if j[0] != i[0]:
                continue
            #set j to be the position
            j = j[1]
            #find the position in the fasta list
            #remember that the fasta sequences start with 1, not 0, so we need to subtract 1 from the position

            try:
                query_base = sequence[j-1]

            except:
                logging.warning(f"position {j} not found in {i[0]}. Skipping...")
                continue

            #extract the context
            context_left = sequence[j-1-context_flank:j-1]
            context_right = sequence[j:j+context_flank]

            context_list.append([context_left, query_base, context_right, i[0]])

    return(context_list)

def drawGenomicContext(context_list: list, show: bool = False, **kwargs: dict) -> None:
    '''
    Draws a sequence logo of a list of genomic contexts.

    Args:
        context_list (list): A list of genomic contexts, where each context is a list of four elements:
            - The sequence of bases to the left of the position.
            - The query base at the position.
            - The sequence of bases to the right of the position.
            - The chromosome name.
        show (bool, optional): Whether to show the plot. Defaults to False.
        **kwargs (dict): Additional keyword arguments to pass to the function. Supported arguments include:
            - save_path (str): The path to save the plot to.

    Returns:
        None

    Raises:
        ImportError: If the required packages matplotlib.pyplot and pandas are not installed.
    '''
    
    import matplotlib.pyplot as plt
    import pandas as pd

    #create a dataframe from the context list
    df = pd.DataFrame(context_list, columns=['context_left', '0', 'context_right', 'chromosome'])

    #split the context into columns by base
    split_left = df['context_left'].apply(lambda x: pd.Series(list(x)))
    split_right = df['context_right'].apply(lambda x: pd.Series(list(x)))

    #rename the columns such that they are numbered. e.g. -3,-2,-1, 0, 1, 2, 3
    split_left.columns = [f"-{i}" for i in range(len(split_left.columns), 0, -1)]
    split_right.columns = [f"{i+1}" for i in range(len(split_right.columns))]

    #make a new dataframe with the context split into columns, query still in the middle
    df_context = pd.concat([split_left, df['0'], split_right], axis=1)

    df_context_freq = df_context.apply(pd.value_counts).fillna(0)              #get the frequencies of each base in each column
    df_context_freq = df_context_freq.div(df_context_freq.sum(axis=0), axis=1) #normalize the frequencies
    df_context_freq = df_context_freq.transpose()
    df_context_freq = df_context_freq*100                                      #set the frequencies to %
    
    #set color scheme for the bases so that they are consistent
    colors = {
                'A':'tab:green', 
                'T':'tab:red',
                'G':'tab:orange',
                'C':'tab:blue',
                'N':'tab:gray',
                }
    
    #plot the dataframe as a stacked barplot
    column_count = len(df_context_freq.columns) #count the number of columns
    df_context_freq.plot.bar(stacked=True, figsize=(column_count+3,4), color = list(colors.values()))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('Frequency (%)')
    plt.xlabel('Relative Position')

    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], bbox_inches='tight', dpi=300)
        logging.info(f"saved plot to {kwargs['save_path']}")

    if show == True:
        plt.show()

    plt.close()
