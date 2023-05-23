# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## validateDNASquence
## imperfectHomologySearch
## findInvertedRepeat
## searchSequenceForRepeats
## saveToJSON

import regex as re

#from .base import validateSquence
from .base import rev_compl


def validateDNASquence(sequence: str) -> bool:
    '''
    Check if the given DNA sequence contains only valid nucleotides.

    Args:
        sequence (str): the DNA sequence to validate.

    Returns:
        bool: True if the sequence contains only valid nucleotides, False otherwise.
    '''

    bases = "ATGCatgc-"
    for i in sequence:
        if i not in bases:
            print("warning, sequence doesn't contain cannonical nucleotides")
            return(False)
    return(True)

def imperfectHomologySearch(sequence: str, query: str, min_homology: float = 0.8, fixed_errors: bool = False, inverted:bool = True) -> list:
    '''
    Search for a query sequence in a given DNA sequence, allowing for a certain number of mismatches/imperfect homology.

    Args:
        sequence (str): the DNA sequence to search in.
        query (str): the query sequence to search for.
        min_homology (float, optional): the minimum homology (similarity) between the query and the found sequence, as a fraction between 0 and 1. Defaults to 0.8.
        fixed_errors (bool, optional): if True, use a fixed number of errors instead of calculating it from min_homology. Defaults to False.
        inverted (bool, optional): if True, search for the reverse complement of the query sequence as well. Defaults to True.

    Returns:
        list: a list of query-match pairs, where each pair is a list containing the query sequence and a list of matching sequences.
    '''
    errors = 0
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

def findInvertedRepeat(sequence: str, query_length: int = 4, min_spacer: int = 4, imperfect_homology: bool = False,
                        min_homology: float = 0.8, fixed_errors: bool = False, inverted: bool = True) -> list:
    '''
    Search for inverted repeats in a given DNA sequence.

    Args:
        sequence (str): the DNA sequence to search in.
        query_length (int, optional): the length of the query sequence to search for. Defaults to 4.
        min_spacer (int, optional): the minimum number of non-matching nucleotides between the two halves of the inverted repeat. Defaults to 4.
        imperfect_homology (bool, optional): if True, allow for a certain number of errors in the inverted repeat. Defaults to False.
        min_homology (float, optional): the minimum homology (similarity) between the query and the found sequence, as a fraction between 0 and 1. Only used if imperfect_homology is True. Defaults to 0.8.
        fixed_errors (bool, optional): if True, use a fixed number of errors instead of calculating it from min_homology. Only used if imperfect_homology is True. Defaults to False.
        inverted (bool, optional): if True, search for the reverse complement of the query sequence as well. Defaults to True.

    Returns:
        list: a list of query-match pairs, where each pair is a list containing the query sequence and a list of matching sequences.
    '''
    query_string=sequence[:query_length]

    if inverted == True:
        query = rev_compl(query_string)
    elif inverted == False:
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
        
def searchSequenceForRepeats(sequence: str, min_query_length: int = 4, max_query_length: int = 25,
                            min_spacer: int = 0, window_size: int = 250, imperfect_homology:bool = False,
                            min_homology:float = 0.8, fixed_errors:bool = False, inverted: bool = True) -> dict:
    '''
    Search a DNA sequence for inverted repeats of various lengths.

    Args:
        sequence (str): the DNA sequence to search in.
        min_query_length (int, optional): the minimum length of the query sequence to search for. Defaults to 4.
        max_query_length (int, optional): the maximum length of the query sequence to search for. Defaults to 25.
        min_spacer (int, optional): the minimum number of non-matching nucleotides between the two halves of the inverted repeat. Defaults to 0.
        window_size (int, optional): the size of the sliding window to use when searching for inverted repeats. Defaults to 250.
        imperfect_homology (bool, optional): if True, allow for a certain number of errors in the inverted repeat. Defaults to False.
        min_homology (float, optional): the minimum homology (similarity) between the query and the found sequence, as a fraction between 0 and 1. Only used if imperfect_homology is True. Defaults to 0.8.
        fixed_errors (bool, optional): if True, use a fixed number of errors instead of calculating it from min_homology. Only used if imperfect_homology is True. Defaults to False.
        inverted (bool, optional): if True, search for the reverse complement of the query sequence as well. Defaults to True.

    Returns:
        dict: a dictionary where the keys are the query sequences and the values are lists of matching sequences.
    '''
    if imperfect_homology:
        print(f"Search has been set to find quasi-pallindromes")
        if fixed_errors:
            print(f"        Allowing up to {fixed_errors} errors/mismatches...")
        if not fixed_errors:
            print(f"        Searching with a minimum of {min_homology} homology")
    if not imperfect_homology:
        print(f"Search has been set to find perfect pallindromes")

    
    list_of_all_pairs=[]
    for query_length in range(min_query_length, max_query_length+1):
        print(f"###Searching for inverted-repeating {query_length}bp-long fragments...")
        sequence_size=len(sequence)
        
        for i in range(sequence_size): #-window_size+1
            seq=sequence[i:i+window_size]
            output_pair_list = findInvertedRepeat(  seq,
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

def saveToJSON(results_dictionary: dict, output_file_name: str = 'json_output.json') -> None:
    '''
    Save a dictionary to a JSON file.

    Args:
        results_dictionary (dict): the dictionary to save.
        output_file_name (str, optional): the name of the output file. Defaults to 'json_output.json'.

    Returns:
        None
    '''
    import json
    with open(output_file_name, 'w') as outfile:
            json.dump(results_dictionary, outfile)

