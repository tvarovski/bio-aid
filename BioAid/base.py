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

def extractSeqFromFastaToList(fasta_file_path):
    '''This function takes a path of a fastafile and extracts all sequence names and
    sequences into a nested list [[title_0, sequence_0], [title_1, sequence_1],...]'''

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

def validateSquence(sequence):
    '''Check if the sequence contains viable nucleotides only'''

    bases = "ATGCatgcN"
    for i in sequence:
        if i not in bases:
            logging.info(f"Sequence doesn't contain cannonical nucleotides: {i}")
            return(False)
        if i == "N":
            logging.warning(f"Warning, sequence contains 'N's")
    return(True)

def compl(base):
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

def rev_compl(seq):
    new_seq = ""
    for base in seq:
        new_base = compl(base)
        new_seq = new_base + new_seq
    return(new_seq)

def unique(list1):
    # this function takes a list1 and returns a list2 with
    # unique elements of list 

    list_set = set(list1) 
    unique_list = (list(list_set))
    return unique_list

def createChrList(chrNum):
    chrList = []
    for i in range(chrNum):
        chrList.append(f"chr{i+1}")
    return(chrList)

def dataFrameImport(directory):
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

def pullGenomicContext(list_of_positions, fasta_file_path, context_flank=5):
    '''This function takes a list of positions [chromosome, position] and a fasta file path and returns a list of
    genomic context of the specified length around each position'''

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

def drawGenomicContext(context_list, show=False, **kwargs):
    '''This function takes a list of genomic contexts and draws a sequence logo of the context.
    If you pass show=True, it will show the plot.
    If you pass save_path=<PATH>, it will save the plot to the path.
    It requires the packages matplotlib.pyplot and panda.s'''
    
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
    df_context_freq.plot.bar(stacked=True, figsize=(column_count+3,4), color=colors.values())
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('Frequency (%)')
    plt.xlabel('Relative Position')

    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], bbox_inches='tight', dpi=300)
        logging.info(f"saved plot to {kwargs['save_path']}")

    if show == True:
        plt.show()

    plt.close()
