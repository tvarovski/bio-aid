# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## findpairings
## findFrequencies
## plot
## runOligoFreqAnalysis

import pandas as pd
from .base import extractSeqFromFastaToList
from .base import validateSquence

def findpairings(sequence: str, pairing_length: int) -> list:
    '''
    Find all oligonucleotides of a given length in a DNA sequence.

    Args:
        sequence (str): the DNA sequence to search in.
        pairing_length (int): the length of the oligonucleotides to find.

    Returns:
        list: a list of all oligonucleotides of length 'pairing_length' found in 'sequence'.
    '''

    pairing_list = []
    for i in range(len(sequence)):
        next_chunk = sequence.upper()[i:i+pairing_length]
        if len(next_chunk) == pairing_length:
            pairing_list.append(next_chunk)
    return(pairing_list)

def findFrequencies(alist: list) -> dict:
    '''
    Count the frequency of each item in a list.

    Args:
        alist (list): the list to count the frequency of.

    Returns:
        dict: a dictionary where the keys are the items in the list and the values are their frequency.
    '''

    sequence_dict = {}

    for i in alist:
        if i in sequence_dict:
            sequence_dict[i]+=1
        else:
            sequence_dict[i] = 1
    return(sequence_dict)

def plot(oligo_dict: dict, title: str) -> None:
    '''
    Plot a histogram of the frequency of each item in a dictionary.

    Args:
        oligo_dict (dict): the dictionary to plot the frequency of.
        title (str): the title of the plot.

    Returns:
        None
    '''
    freq_list = list(oligo_dict.items())
    freq_list.sort(key=lambda x:x[1], reverse=True)
    df = pd.DataFrame(freq_list, columns=['oligo', 'frequency'])
    total_oligos = df.frequency.sum()
    df['frequency'] = df.frequency / total_oligos
    df.plot(kind='bar', x='oligo', figsize=(15,8), title=title+" Oligonucleotide Frequency", xlabel='Oligonucleotide', ylabel='Frequency (%)')

def runOligoFreqAnalysis(file_path: str) -> None:
    '''
    Run oligonucleotide frequency analysis on all sequences in a FASTA file.

    Args:
        file_path (str): the path to the FASTA file to analyze.

    Returns:
        None
    '''
    import matplotlib.pyplot as plt
    fa_seq_list = extractSeqFromFastaToList(file_path)

    for fa_seq in fa_seq_list:
        sequence = fa_seq[1]
        title = fa_seq[0]

        if not validateSquence(sequence):
            print(f'skipping "{title}", sequence validation failed')
            continue

        pairing_list_1 = findpairings(sequence,2)
        seqdeck_1 = findFrequencies(pairing_list_1)

        pairing_list_2 = findpairings(sequence,1)
        seqdeck_2 = findFrequencies(pairing_list_2)

        plot(seqdeck_1, title)
        plot(seqdeck_2, title)
        plt.show()


if __name__ == "__main__":
    import sys
    file_path = sys.argv[1]
    runOligoFreqAnalysis(file_path)

