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

def wordsInSequence(sequence, treeLevel):
  # takes a genetic sequence and a treeLevel (word length) and returns
  # a list with words in that sequence

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

def movingWindow(sequence, treeLevel=8, chunkSize=20, complexity_threshold=0.2, showGraph=False):
  # takes a genetic sequence and calculates the linguistic complexity scores
  # for windows of size chunkSize for word lengths from 1 to treeLevel.
  # returns the lowest score window in format:
  # [Boolean, [complexity, treeLevel, sequence]] where Boolean denotes if the
  # complexity score is lower (True) than complexity_threshold
  # optional: draw the complexity graph for the length of the sequence
  # by setting showGraph = True

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
    ax.set_title(lowest[2])

  if lowest[0] < complexity_threshold:
    return([True, lowest[0], lowest[1]])
    #return([True, lowest])
  else:
    return([False, lowest[0], lowest[1]])
    #return([False, lowest])

def createLogFile(seqlist, filename, complexity_threshold = 0.2, chunkSize=20):
  # takes a list of sequences for analysis, the name of the output file,
  # complexity threshold, and windowsize (chunkSize), and outputs a log file
  # with lowesst complexity score for each sequence
  outputf = open(filename, "w")
  truemmbir = 0
  for i in seqlist:
    seqScore = movingWindow(i, complexity_threshold = complexity_threshold, chunkSize = chunkSize)
    if seqScore[0] == False:
      truemmbir += 1
    outputf.write(str(seqScore[1])+"\n")
  print("sequences above threshold:", truemmbir)

  outputf.close()

def createSequenceList(filename):
  # takes a filename and extracts sequences from it, and outputs them to a list
  # has some checking (all sequences shorter than 100bp will be discarted)
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

def createDataFrameFromLogFile(logfile):
  # creates a pandas dataframe out of the complexity scores log file
  logfile = open(logfile, "r")
  lines = logfile.readlines()
  loglist=[]
  for line in lines:
    if "[" in line:
        loglist.append(line.strip("\n"))

  loglist2=[]
  for i in range(len(loglist)):
    loglist[i] = ast.literal_eval(loglist[i]) 
  print("total sequences:", len(loglist))

  df = DataFrame(loglist,columns=['complexity', 'level', 'sequence'])
  
  return df