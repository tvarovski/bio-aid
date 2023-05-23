# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import regex as re

def createMMBSearchReference(file_path_in: str, file_path_out: str) -> None:
    """
    Reads a FASTA file containing genome data and creates a new file with a modified header line for each chromosome.
    It's purpose is to create a FASTA file that can be used as a reference for MMBSearch.

    Args:
        file_path_in (str): The path to the input FASTA file.
        file_path_out (str): The path to the output file to be created.

    Returns:
        None
    """
    #good for human genome
    file_in = open(file_path_in, "r")
    save=False
    readlines=True
    while readlines:
        try:
            line = file_in.readline()
        except:
            print("Failed to read line. EOF? Exiting...")
            readlines=False
            break

        if line[0] == ">":
            if re.search("chromosome.*Primary.Assembly$", line) != None:
                print(line)
                linewords = line.split()
                chromosome=linewords[4].strip(",")
                if (chromosome=="1") | (chromosome=="2"):
                    chromosome=f">chr0{chromosome}\n"
                elif chromosome=="X":
                    chromosome=f">chrX\n"
                elif chromosome=="Y":
                    chromosome=f">chrY\n"
                else:
                    chromosome=f">chr{chromosome}\n"
                line=chromosome
                print(line)
                save=True
            #fix for mitochondrial chr
            elif re.search("mitochondrion, complete genome$", line) != None:
                print(line)
                line=f">chrM\n"
                print(line)
                save=True
            else:
                save=False

        if save:
            with open(file_path_out, "a") as file_out:
                file_out.write(line)