# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## verifyImperfectHomology
## add_gene
## add_exon
## createAnnotatedOutput

import pandas as pd
from pyensembl import EnsemblRelease
from pyensembl.exon import Exon
import logging
import regex as re
from ..base import createChrList
from ..base import rev_compl
from ..complexity import movingWindow


def verifyImperfectHomology(ref, query, min_homology=0.8):

  mmbir = rev_compl(query)
  mmbir_errors = round(len(mmbir)*(1-min_homology))
  if mmbir_errors > 10:
    mmbir_errors = 10
  output_list = re.findall( '(' + mmbir + '){e<=' + str(mmbir_errors) + '}', ref)
  print(f'Found {len(output_list)} possible templates for {query}: {output_list}')

  if len(output_list) > 0:
    return(True)
  else:
    return(False)

def add_gene(row):

  data = EnsemblRelease(104)
  #print("loaded genome! #addgene")
  output = data.gene_names_at_locus(contig=row['chr'].strip("chr"), position=int(row["iBirStart"]))
  print(output)
  #output = str(output)
  out_list = []
  for gname in output:
    if gname != "":
      out_list.append(gname)
      #print("appended_gene!")

  return(out_list) #output

def add_exon(row):

  data = EnsemblRelease(104)
  #print("loaded genome! #addexon")
  output = data.exons_at_locus(contig=row['chr'].strip("chr"), position=int(row["iBirStart"]))
  #output = str(output)
  out_list = []
  for exon in output:
    if isinstance(exon, Exon):
      out_list.append(exon.exon_id)

  return(out_list) #output

def createAnnotatedOutput(f_name, path, output_f_name):
  import os

  chrList=createChrList(25)
  df = pd.DataFrame()

  for chr in chrList:
    try:
      if os.stat(f"{path}{chr}/{f_name}").st_size == 0:
        print(f"WARNING! File {path}{chr}/{f_name} is empty. Skipping...")
        continue
      with open(f"{path}{chr}/{f_name}") as f:
        lines = f.readlines()
    except:
      print(f"WARNING! Couldn't read {path}{chr}/{f_name}. Make sure it's there.")
      lines = None
      continue

    if (type(lines) == None):
      print("WARNING! the type of lines is None. Skipping...")

    else:
      print(f"hello {chr}")
      for i in range(len(lines)):
        if "###" in lines[i]:
          i-=1
          if "Microhomology Insertion:" in lines[i+11]:
            var0 = lines[i].strip()
            var1 = lines[i+1].strip()
            var2 = lines[i+2].strip("iBirStart:").strip()
            var3 = lines[i+3].strip("Consensus/cluster number:").strip()
            var4 = lines[i+4].strip("iDepth:").strip()
            var5 = lines[i+5].strip("sBir:").strip()
            var6 = lines[i+6].strip("sBirReversed:").strip()
            var7 = lines[i+7].strip()
            var8 = lines[i+8].strip("ref:").strip()
            var9 = lines[i+9].strip("bir:").strip()
            var10 = lines[i+10].strip("TEM:").strip()
            var11 = lines[i+11].strip("Microhomology Insertion:").strip()
            var12 = lines[i+12].strip("Microhomology template:").strip()
            var13 = lines[i+13].strip("readStart: ").strip()
            var14 = lines[i+14].strip("ref:").strip()
            var15 = lines[i+15].strip("bir:").strip()
            var16 = lines[i+16].strip()
            var17 = lines[i+17].strip()
            var18 = lines[i+18].strip()
            var19 = lines[i+19].strip()
          else:
            var0 = lines[i].strip()
            var1 = lines[i+1].strip()
            var2 = lines[i+2].strip("iBirStart:").strip()
            var3 = lines[i+3].strip("Consensus/cluster number:").strip()
            var4 = lines[i+4].strip("iDepth:").strip()
            var5 = lines[i+5].strip("sBir:").strip()
            var6 = lines[i+6].strip("sBirReversed:").strip()
            var7 = lines[i+7].strip()
            var8 = lines[i+8].strip("ref:").strip()
            var9 = lines[i+9].strip("bir:").strip()
            var10 = lines[i+10].strip("TEM:").strip()
            var11 = "Empty"
            var12 = "Empty"
            var13 = lines[i+11].strip("readStart: ").strip()
            var14 = lines[i+12].strip("ref:").strip()
            var15 = lines[i+13].strip("bir:").strip()
            var16 = lines[i+14].strip()
            var17 = lines[i+15].strip()
            var18 = lines[i+16].strip()

          event_dict = {
                        "chr":chr,
                        "iBirStart":var2,
                        "Consensus/cluster_number":var3,
                        "iDepth":var4,
                        "sBir":var5,
                        "sBirReversed":var6,
                        "ref":var8,
                        "bir":var9,
                        "TEM":var10,
                        "Microhomology_Insertion":var11,
                        "Microhomology_template":var12,
                        "readStart":var13,
                        "var16":var16,
                        "var17":var17,
                        "var18":var18
                        }
          
          #make sure that the ref and bir are not empty, otherwise the movingWindow will fail
          if (len(var8) == 0) or (len(var9) == 0):
            print(f"ref or bir are empty for {path}{chr}/{f_name}")
            print("Cheeck if they run correctly in the original output file...")
            print("Skipping this event...")

            #print the same messages to the error log
            logging.basicConfig(filename='error.log', level=logging.DEBUG)
            logging.debug(f"ref or bir are empty for {path}{chr}/{f_name}")
            logging.debug("Cheeck if they run correctly in the original output file...")
            logging.debug("Skipping this event...")
            continue

          event_df = pd.DataFrame.from_records([event_dict])
          df = pd.concat([df, event_df], ignore_index=True)

  print("Starting complexity")

  try:
    df[['ref_complexity_fail', 'ref_complexity_score', 'ref_complexity_tree']] = df.apply(lambda row: movingWindow(row.ref, complexity_threshold=0.2), axis=1, result_type ='expand')
  except:
    print(f"Error with ref_complexity for {path}{chr}/{f_name}")

  try:
    df[['bir_complexity_fail', 'bir_complexity_score', 'bir_complexity_tree']] = df.apply(lambda row: movingWindow(row.bir, complexity_threshold=0.2), axis=1, result_type ='expand')
  except:
    print(f"Error with bir_complexity for {path}{chr}/{f_name}")

  print("Starting homology check")
  df['homology_check_ref'] = df.apply(lambda row: verifyImperfectHomology(row.ref, row.sBir), axis=1)
  df['homology_check_bir'] = df.apply(lambda row: verifyImperfectHomology(row.bir, row.sBir), axis=1)
  print("Starting gene/exon annotation")
  df["genes"] = df.apply(lambda row: add_gene(row), axis=1)
  df["exones"] = df.apply(lambda row: add_exon(row), axis=1)

  df.to_csv(output_f_name, sep="\t",index=False)
  print(f"finished annotation... DF saved to {output_f_name}")