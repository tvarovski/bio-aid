# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import regex as re

def createMMBSearchReference(file_path_in, file_path_out):
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