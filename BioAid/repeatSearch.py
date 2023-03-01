# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## validateDNASquence
## imperfectHomologySearch
## findInvertedRepeat
## searchSequenceForRepeats
## saveToJSON

import regex as re
import json

#from .base import validateSquence
from .base import rev_compl


def validateDNASquence(sequence):
  '''Check if the sequence contains viable nucleotides only'''

  bases = "ATGCatgc-"
  for i in sequence:
    if i not in bases:
      print("warning, sequence doesn't contain cannonical nucleotides")
      return(False)
  return(True)

def imperfectHomologySearch(sequence, query, min_homology=0.8, fixed_errors=False, inverted=True):

  if min_homology:
    errors = round(len(query)*(1-min_homology))
  if fixed_errors:
    errors = fixed_errors
  #print(f"searching with {errors} errors...")

  output_list = re.findall( '(' + query + '){e<=' + str(errors) + '}', sequence)

  if inverted == True:
    query=rev_compl(query)
  
  query_match_pairs=[query, output_list]
  if len(output_list) > 0:
    if len(output_list) > 1:
      print("This is unusual...")
    print(f'Found possible template(s) for {query}: {output_list}')
    return(query_match_pairs)
  else:
    return([])

def findInvertedRepeat(sequence, 
                       query_length=4,
                       min_spacer=4,
                       imperfect_homology=False,
                       min_homology=0.8,
                       fixed_errors=False,
                       inverted=True):

  query_string=sequence[:query_length]

  if inverted == True:
    query = rev_compl(query_string)
  if inverted == False:
    query = query_string
    
  sequence=sequence[query_length+min_spacer:]
  query_complement = rev_compl(query)
  output_pair_list=[]

  for i in range(len(sequence)-query_length+1):
    #print(sequence[i:query_length+i])
    if imperfect_homology:
      
      
      output_pair=imperfectHomologySearch(sequence[i:query_length+i],
                                          query,
                                          min_homology=min_homology,
                                          fixed_errors=fixed_errors,
                                          inverted=inverted)

      if output_pair != []:
        output_pair_list.append(output_pair)

    else:
      if sequence[i:query_length+i] == query:

        if inverted == True:
          print(f"Success {query_complement}, {sequence[i:query_length+i]}")
          output_pair_list.append([query_complement, [sequence[i:query_length+i]]])

        if inverted == False:
          print(f"Success {query}, {sequence[i:query_length+i]}")
          output_pair_list.append([query, [sequence[i:query_length+i]]])

  return(output_pair_list)
    
def searchSequenceForRepeats(sequence,
                             min_query_length=4,
                             max_query_length=25,
                             min_spacer=0,
                             window_size=250,
                             imperfect_homology=False,
                             min_homology=0.8,
                             fixed_errors=False,
                             inverted=True):
  
  if imperfect_homology:
    print(f"Search has been set to find quasi-pallindromes")
    if fixed_errors:
      print(f"    Allowing up to {fixed_errors} errors/mismatches...")
    if not fixed_errors:
      print(f"    Searching with a minimum of {min_homology} homology")
  if not imperfect_homology:
    print(f"Search has been set to find perfect pallindromes")

  
  list_of_all_pairs=[]
  for query_length in range(min_query_length, max_query_length+1):
    print(f"###Searching for inverted-repeating {query_length}bp-long fragments...")
    sequence_size=len(sequence)
    
    for i in range(sequence_size): #-window_size+1
      seq=sequence[i:i+window_size]
      output_pair_list = findInvertedRepeat(seq,
                                            query_length=query_length,
                                            min_spacer=min_spacer,
                                            imperfect_homology=imperfect_homology,
                                            min_homology=min_homology,
                                            fixed_errors=fixed_errors,
                                            inverted=inverted)
      for pair in output_pair_list:
        list_of_all_pairs.append(pair)

  results_dictionary={}
  print("Consolodating Results...")
  for pair in list_of_all_pairs:
    query=pair[0]
    if query in results_dictionary:
      if pair[1][0] not in results_dictionary[query]:
        results_dictionary[query].append(pair[1][0])
    else:
      results_dictionary[query] = [pair[1][0]]
  print("Results were consolidated...")
  return(results_dictionary)

def saveToJSON(results_dictionary, output_file_name='json_output.json'):

  with open(output_file_name, 'w') as outfile:
      json.dump(results_dictionary, outfile)

