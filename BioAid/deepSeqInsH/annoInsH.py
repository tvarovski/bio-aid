# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## dataFrameImport
## matchSequence
## findReadLength
## trimseq
## classify

import regex as re
import os
import pandas as pd

def dataFrameImport(directory,datatype = 'sam'):
  import pysam

  #This function imports all files inside of the specified DIR path into two lists; of panda DF and sample names
  frames_list_samples = []
  sample_names = []

  counter=1

  #import all files in the directory
  for filename in os.listdir(directory):
      if filename.endswith(f".{datatype}"):
          sample_names.append(filename[:8].strip())
          path = os.path.join(directory, filename)
          globals()[f"sample{counter}"] = pd.DataFrame(({'name': x.qname, 'seq': x.seq, 'qual': x.query_qualities} for x in pysam.Samfile(path).fetch()))
          print(f'dataframe {sample_names[-1]} has {len(globals()[f"sample{counter}"])} total rows')
          frames_list_samples.append(globals()[f"sample{counter}"])
          counter+=1
  print(f"found {len(frames_list_samples)} samples in {directory}")
  print(len(sample_names), sample_names)

  return(frames_list_samples, sample_names)

def matchSequence(row, sequence):
  #for searching for perfectly matched sequences.
  if (sequence in row):
    return(True)
  else:
    return(False)

def findReadLength(row):
  #returns a length of a read
  return(len(row))

def trimseq(row, f_seq, r_seq, consensus):
  #trims sequences to the 5' end of primers and finds reads that support non-deletions.
  #allows for 2 errors in non-deletions.

  if row.forward_primer == True:
    f_seq_index = row.seq.find(f_seq)
    if f_seq_index == -1:
      print("no forward primer found")
      row.forward_primer = False
    row.seq = row.seq[f_seq_index:]

  if row.reverse_primer == True:
    r_seq_index = row.seq.find(r_seq)
    if r_seq_index == -1:
      print("no reverse primer found")
      row.reverse_primer = False
    row.seq = row.seq[:r_seq_index+len(r_seq)]
  if (row.forward_primer == False) & (row.reverse_primer == False):
    row.seq = "no_primers"

  if re.search('('+row.seq+')'+'{e<=2}', consensus):
    row.seq = "no_excision"
  return(row)

def classify(row, consensus, classification):
  #this function classifies the sequence based on the provided consensus, allows for 2 errors.

  if row.classification != "other":
    return row

  if re.search('('+row.seq+')'+'{e<=2}', consensus):
  #if row.seq in consensus:
    row.classification = classification
  return row