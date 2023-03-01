# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## findpairings
## findFrequencies
## plot
## runOligoFreqAnalysis

import matplotlib.pyplot as plt
import pandas as pd
from .base import extractSeqFromFastaToList
from .base import validateSquence

def findpairings(sequence, pairing_length):
  '''finds all oligonucleotides in 'sequence' of length 'pairing_length' based on 
  a sliding window and returns a list that contains them'''

  pairing_list = []
  for i in range(len(sequence)):
    next_chunk = sequence.upper()[i:i+pairing_length]
    if len(next_chunk) == pairing_length:
      pairing_list.append(next_chunk)
  return(pairing_list)

def findFrequencies(alist):
  '''takes a list of oligonucleotides and returns a dictionary where 
  oligonucleotides are keys and their count in the list are values'''

  sequence_dict = {}

  for i in alist:
    if i in sequence_dict:
      sequence_dict[i]+=1
    else:
      sequence_dict[i] = 1
  return(sequence_dict)

def plot(oligo_dict, title):
  '''this function plots the histogram of frequencies based on the dict output of 
  'findFrequencies' function'''

  freq_list = list(oligo_dict.items())
  freq_list.sort(key=lambda x:x[1], reverse=True)
  df = pd.DataFrame(freq_list, columns=['oligo', 'frequency'])
  total_oligos = df.frequency.sum()
  df['frequency'] = df.frequency / total_oligos
  df.plot(kind='bar', x='oligo', figsize=(15,8), title=title+" Oligonucleotide Frequency", xlabel='Oligonucleotide', ylabel='Frequency (%)')

def runOligoFreqAnalysis(file_path):
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



#file_path = sys.argv[1]
#runOligoFreqAnalysis(file_path)
