

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
  print(f"Extraction of sequence information from {fasta_file_path} finished.")
  return fasta_list


def validateSquence(sequence):
  '''Check if the sequence contains viable nucleotides only'''

  bases = "ATGCatgc"
  for i in sequence:
    if i not in bases:
      print("warning, sequence doesn't contain cannonical nucleotides")
      return(False)
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