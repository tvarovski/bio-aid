# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## wordsInSequence
## findComplexity
## movingWindow
## createLogFile
## createSequenceList
## createDataFrameFromLogFile

import pandas as pd
import ast
from pandas import DataFrame
from .base import unique

def wordsInSequence(sequence: str, treeLevel: int) -> list:
    """
    Extracts all possible words of a specified length from a genetic sequence.

    Args:
        sequence (str): The genetic sequence to extract words from.
        treeLevel (int): The length of the words to extract.

    Returns:
        list: A list of words, where each word is a string of length treeLevel.

    Examples:
        >>> wordsInSequence("ATCGATCG", 3)
        ['ATC', 'TCG', 'CGA', 'GAT', 'ATC', 'TCG']

        >>> wordsInSequence("ATCGATCG", 4)
        ['ATCG', 'TCGA', 'CGAT', 'GATC', 'ATCG']
    """
    length = len(sequence)
    max_possible_words_in_sequence = length-(treeLevel-1)
    wordList = []
    for i in range(max_possible_words_in_sequence):
        wordList.append(sequence[0+i:treeLevel+i])
    return wordList

def findComplexity(sequence, treeLevel, complexity_threshold):
    # takes a genetic sequence and calculates the linguistic complexity 
    # for that sequence for word length treeLevel

    wordList = wordsInSequence(sequence,treeLevel)
    wordList = unique(wordList)

    if len(sequence) > 4**treeLevel:
        complexity = len(wordList)/(4**treeLevel)
    else:
        complexity = len(wordList)/(len(sequence)-(treeLevel-1))

    if (complexity < complexity_threshold) & (complexity > 0.2):
        print("Complexity at tree level " + str(treeLevel) + " is " + str(complexity) + " for sequence: "+str(sequence))

    return ([complexity, treeLevel, sequence])

def movingWindow(sequence: str, treeLevel: int = 8, chunkSize: int = 20, complexity_threshold: float = 0.2, showGraph: bool = False) -> list:
    """
    Calculates the linguistic complexity scores for windows of a specified size and word lengths from 1 to a specified level.

    Args:
        sequence (str): The genetic sequence to calculate complexity scores for.
        treeLevel (int, optional): The maximum length of words to use in the complexity calculation. Defaults to 8.
        chunkSize (int, optional): The size of the windows to use in the complexity calculation. Defaults to 20.
        complexity_threshold (float, optional): The threshold below which a window is considered to have low complexity. Defaults to 0.2.
        showGraph (bool, optional): Whether to show a graph of the complexity scores. Defaults to False.

    Returns:
        list: A list containing the lowest complexity score window in the format [Boolean, complexity, treeLevel], where Boolean denotes if the complexity score is lower (True) than complexity_threshold.

    Raises:
        ValueError: If the sequence is too short to calculate complexity scores.

    Examples:
        >>> movingWindow("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", treeLevel=4, chunkSize=10, complexity_threshold=0.1, showGraph=True)
        [True, 0.0, 4]

        >>> movingWindow("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", treeLevel=10, chunkSize=5, complexity_threshold=0.5, showGraph=True)
        [False, 0.5, 10]
    """

    complexityList=[]

    for i in range(1, treeLevel+1):
        for chunk in range(len(sequence)-(chunkSize-1)):

            chunkSeq = sequence[chunk:(chunk+chunkSize)]
            complexityList.append(findComplexity(chunkSeq, i, complexity_threshold))

    try:
        lowest=min(complexityList, key=lambda x: x[0])

    except ValueError:
        print(f"ERROR with sequence: {sequence}. It is probably too short len={len(sequence)}. Skipping it.")
        print("complexityList: ", complexityList)
        lowest = [0,0,0]

    #print(complexityList)

    lowLvl=lowest[1]

    if showGraph:

        df = DataFrame(complexityList,columns=['complexity','tree','sequence'])
        df = df[['complexity','tree']]
        ax = df[df.tree == lowLvl].filter(items=["complexity"]).plot(fontsize=15, grid=True, figsize=(20,8))
        ax.axhline(complexity_threshold, color="red", linestyle="--")
        ax.axhline(1, color="green", linestyle="--")
        ax.set_ylim(ymin=0)
        ax.set_title(str(lowest[2]))

    if lowest[0] < complexity_threshold:
        return([True, lowest[0], lowest[1]])
        #return([True, lowest])
    else:
        return([False, lowest[0], lowest[1]])
        #return([False, lowest])

def createLogFile(seqlist: list, filename: str, complexity_threshold: float = 0.2, chunkSize: int = 20) -> None:
    """
    Creates a log file containing the lowest linguistic complexity score for each sequence in a list.

    Args:
        seqlist (list): A list of genetic sequences to calculate complexity scores for.
        filename (str): The name of the output file to write the log to.
        complexity_threshold (float, optional): The threshold below which a sequence is considered to have low complexity. Defaults to 0.2.
        chunkSize (int, optional): The size of the windows to use in the complexity calculation. Defaults to 20.

    Returns:
        None
    """
    outputf = open(filename, "w")
    truemmbir = 0
    for i in seqlist:
        seqScore = movingWindow(i, complexity_threshold = complexity_threshold, chunkSize = chunkSize)
        if seqScore[0] == False:
            truemmbir += 1
        outputf.write(str(seqScore[1])+"\n")
    print("sequences above threshold:", truemmbir)

    outputf.close()

def createSequenceList(filename: str) -> list:
    """
    Extracts genetic sequences from a file and returns them as a list.

    Args:
        filename (str): The name of the file to extract sequences from.

    Returns:
        list: A list of genetic sequences, where each sequence is a string.
    """
    seqfile = open(filename, "r")
    lines = seqfile.readlines()
    seqlist=[]
    for line in lines:
        if "bir:" in line:
            if len(line) > 100:
                seqlist.append(line.strip("\n").strip("bir:-"))
    seqlist = list(dict.fromkeys(seqlist))
    print(len(seqlist))

    longseq=0
    for i in seqlist:
        if len(i) > 300:
            longseq+=1
    print(longseq)

    return seqlist

def createDataFrameFromLogFile(logfile_path: str) -> DataFrame:
    """
    Creates a pandas DataFrame from a log file containing linguistic complexity scores.

    Args:
        logfile (str): The name of the log file to create the DataFrame from.

    Returns:
        DataFrame: A pandas DataFrame containing the complexity scores, where each row represents a sequence and contains the complexity score, the word length, and the sequence itself.
    """
    logfile = open(logfile_path, "r")
    lines = logfile.readlines()
    loglist=[]
    for line in lines:
        if "[" in line:
                loglist.append(line.strip("\n"))

    for i in range(len(loglist)):
        loglist[i] = ast.literal_eval(loglist[i]) 
    print("total sequences:", len(loglist))

    df = DataFrame(loglist,columns=['complexity', 'level', 'sequence'])
    
    return df