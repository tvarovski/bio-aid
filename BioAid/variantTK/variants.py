# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

## Contents:
## qualitySNPFilter
## filterByAD
## findDominantAF
## findType
## findZygosity
## filterFromClones
## drawSNPMap
## findGenomeLength
## generateRandomLoci
## findChrInLinearGenome
## findSpectraStrandwise
## renameChrToRoman
## draw_random_SNPMap
## drawCombinedSNPMap

import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from natsort import index_natsorted
import seaborn as sns
from numpy.core.numeric import NaN


def qualitySNPFilter(frames_list):

  for i in range(len(frames_list)):
    frames_list[i] = frames_list[i][frames_list[i].TYPE =='SNP']
    frames_list[i] = frames_list[i].rename(columns=lambda c: 'AF' if 'AF' in c else c)
    frames_list[i] = frames_list[i].rename(columns=lambda c: 'AD' if 'AD' in c else c)
    
    #this line may be problematic for a good result and if not needed comment out
    frames_list[i]['AF'] = frames_list[i].apply(lambda row: findDominantAF(row), axis=1)

    frames_list[i] = frames_list[i].astype({'AF': float})
    frames_list[i] = frames_list[i][frames_list[i].AF >= 0.35]
    #also filter out 1.5n variants?

    frames_list[i]['EVAL_AD'] = frames_list[i].apply(lambda row: filterByAD(row), axis=1)
    frames_list[i] = frames_list[i][frames_list[i].EVAL_AD == 'True']
    frames_list[i]['SPECTRA'] = frames_list[i].apply(lambda row: findType(row.REF, row.ALT), axis=1)
  return frames_list

def filterByAD(row):
  reads=row['AD'].split(",")
  if int(reads[1]) < 5:
    return "False"
  if int(reads[0])+int(reads[1]) > 9:
    return "True"
  else:
    return "False"

def findDominantAF(row):
  reads=row['AF'].split(",")
  return(max(reads))

def findType(colA, colB):
  if ((colA == 'C') & (colB == 'T') | (colA == 'G') & (colB == 'A')):
    return('C_to_T')
  if ((colA == 'C') & (colB == 'A') | (colA == 'G') & (colB == 'T')):
    return('C_to_A')
  if ((colA == 'C') & (colB == 'G') | (colA == 'G') & (colB == 'C')):
    return('C_to_G')
  if ((colA == 'T') & (colB == 'C') | (colA == 'A') & (colB == 'G')):
    return('T_to_C')
  if ((colA == 'T') & (colB == 'G') | (colA == 'A') & (colB == 'C')):
    return('T_to_G')
  if ((colA == 'T') & (colB == 'A') | (colA == 'A') & (colB == 'T')):
    return('T_to_A')

def findZygosity(row):
  if (row["AF"] >= 0.85):
    return('Homozygous')
  else:
    return('Heterozygous')

def filterFromClones(frames_list):
  #this will subtract all common rows between normal/control and experimental
  #sample. I might need to revist to filter variants by position only not all columns
  subtracted_df_list = []
  for df in frames_list:
    temp_df = df.drop(columns=['AD','AF'])
    len_original=len(temp_df)
    for df_2 in frames_list:
      if not df.equals(df_2):
        df_2=df_2.drop(columns=['AD','AF'])
        temp_df = pd.merge(temp_df,df_2, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)

    len_final=len(temp_df)
    print(f"removed {len_original-len_final} rows from df")
    subtracted_df_list.append(temp_df)
  return subtracted_df_list

def drawSNPMap(pd_df, df_chr_lengths, chr_starts_df, title, sample_names, saveMap=True):

  pd_df = pd_df.append(chr_starts_df, ignore_index=True)

  fig, ax = plt.subplots()

  df_chr_lengths.plot(kind='barh',legend=False, ax=ax, x="chromosome", y="end_position", fontsize=15, figsize=(20,6), edgecolor="black", linewidth=1, color="beige", zorder=0, title=title, label="Position")
  
  pd_df = pd_df.sort_values(by="CHROM", key=lambda x: np.argsort(index_natsorted(pd_df["CHROM"])))
  
  #To rename legend elements, change plot markers/colors, modify here
  label_dict = { "G": "G->N", "C": "C->N", "A": "A->N", "T": "T->N", "REF": "Reference Allele"}
  markers = {"Homozygous": "x", "Heterozygous": "|"}
  palette = {"G": "green", "C": "red", "A": "gray", "T": "gray"}

  sns.scatterplot(ax=ax,data=pd_df, x="POS", y="CHROM", hue="REF", palette=palette, style="ZYGOSITY", markers=markers, s=120, alpha=0.7,zorder=1, linewidth=1.5)
  
  ax.set_xlabel("POSITION")
  ax.set_ylabel("CHROMOSOME")

  legend_texts = plt.legend().get_texts()
  for i in range(len(legend_texts)):
    label_str = legend_texts[i].get_text()
    if label_str in label_dict:
      new_label = label_dict.get(label_str)

      legend_texts[i].set_text(new_label)

  if saveMap:
    plt.savefig('figures/'+title+'.png', transparent=False)

def findGenomeLength(chromosomes):

  genome_length = 0
  for chr in chromosomes:
    genome_length+=chr[1]

  return(genome_length)

def generateRandomLoci(n, chromosomes):

  genome_length=findGenomeLength(chromosomes)
  return(random.sample(range(0, genome_length), n))

def findChrInLinearGenome(chromosomes, random_loci):

  output_list=[]
  for locus in random_loci:
    linear_genome_slider=0
    linear_genome_slider_previous=0
    for chromosome in chromosomes:
      linear_genome_slider+=chromosome[1]

      if locus < linear_genome_slider:
        #print(f'slider:{linear_genome_slider}, chr:{chromosome[0]}, locus:{locus}')
        output_list.append([chromosome[0], locus-linear_genome_slider_previous])
        linear_genome_slider_previous+=chromosome[1]
        break
      linear_genome_slider_previous+=chromosome[1]

  return(output_list)

def findSpectraStrandwise(colA, colB):
  if ((colA == 'C') & (colB == 'T')):
    return('C_to_T')
  elif ((colA == 'C') & (colB == 'A')):
    return('C_to_A')
  elif ((colA == 'C') & (colB == 'G')):
    return('C_to_G')
  elif ((colA == 'T') & (colB == 'C')):
    return('T_to_C')
  elif ((colA == 'T') & (colB == 'G')):
    return('T_to_G')
  elif ((colA == 'T') & (colB == 'A')):
    return('T_to_A')

  elif ((colA == 'G') & (colB == 'A')):
    return('G_to_A')
  elif ((colA == 'G') & (colB == 'T')):
    return('G_to_T')
  elif ((colA == 'G') & (colB == 'C')):
    return('G_to_C')
  elif ((colA == 'A') & (colB == 'G')):
    return('A_to_G')
  elif ((colA == 'A') & (colB == 'C')):
    return('A_to_C')
  elif ((colA == 'A') & (colB == 'T')):
    return('A_to_T')

def renameChrToRoman(df):
  roman_names={"chr1": "I",
               "chr2": "II",
               "chr3": "III",
               "chr4": "IV",
               "chr5": "V",
               "chr6": "VI",
               "chr7": "VII",
               "chr8": "VIII",
               "chr9": "IX",
               "chr10": "X",
               "chr11": "XI",
               "chr12": "XII",
               "chr13": "XIII",
               "chr14": "XIV",
               "chr15": "XV",
               "chr16": "XVI"}
  df.replace({"chromosome": roman_names}, inplace=True)
  return(df)

def draw_random_SNPMap(pd_df, df_chr_lengths, chr_starts_df, title, saveMap=True):

  #To rename legend elements, change plot markers/colors, modify here

  label_dict = {"cen": 'Centromere'}
  markers = {"random": "$|$", "cen": "o"}
  palette = {"random": "black", "cen": "white"}

  
  pd_df = pd_df.append(chr_starts_df, ignore_index=True)
  pd_df = pd_df.sort_values(by="CHROM", key=lambda x: np.argsort(index_natsorted(pd_df["CHROM"])))
  renameChrToRoman(pd_df)


  fig, ax = plt.subplots()

  renameChrToRoman(df_chr_lengths)

  df_chr_lengths.plot(kind='barh',legend=False, ax=ax, x="chromosome", y="end_position", fontsize=40, figsize=(80,24), edgecolor="black", linewidth=1, color="lightgray", zorder=0, label="Position")
  sns.scatterplot(ax=ax, data=pd_df, x="POS", y="CHROM", style="REF", markers=markers, hue="REF", palette=palette, s=1500, alpha=1,zorder=1, linewidth=1, edgecolor="black", legend=False)

  ax.set_xlabel("POSITION", fontsize=40)
  ax.set_ylabel("CHROMOSOME", fontsize=40)
  fig.suptitle(title, fontsize=60)

  legend_texts = plt.legend().get_texts()
  for i in range(len(legend_texts)):
    label_str = legend_texts[i].get_text()
    if label_str in label_dict:
      new_label = label_dict.get(label_str)

      legend_texts[i].set_text(new_label)

  if saveMap:
    plt.savefig('figures/'+title+'.png', transparent=False, dpi=200)

def drawCombinedSNPMap(pd_df, df_chr_lengths, chr_starts_df, title, sample_names, saveMap=True):

  #To rename legend elements, change plot markers/colors, modify here
  label_dict = { "G": "G->N",
                 "C": "C->N",
                 "A": "A->N",
                 "T": "T->N",
                 "SPECTRA_STRANDWISE": "Reference Allele",
                 "cen": 'Centromere'}
                
  markers = {"C_to_T": "$|$",
             "C_to_A": "$|$",
             "C_to_G": "$|$",
             "G_to_A": "$|$",
             "G_to_C": "$|$",
             "G_to_T": "$|$",
             'T_to_A': "$|$",
             'T_to_C': "$|$",
             'A_to_C': "$|$",
             'A_to_G': "$|$",
             'A_to_T': "$|$",
             'T_to_G': "$|$",
             'cen': "o"}

  palette = {"C_to_T": "black",
             "C_to_A": "black",
             "C_to_G": "black",
             "G_to_A": "red",
             "G_to_C": "blue",
             "G_to_T": "olive",
             'T_to_A': "gray",
             'T_to_C': "gray",
             'A_to_C': "gray",
             'A_to_G': "gray",
             'A_to_T': "gray",
             'T_to_G': "gray",
             'cen': "white"}

  pd_df = pd_df.append(chr_starts_df, ignore_index=True)
  pd_df = pd_df.sort_values(by="CHROM", key=lambda x: np.argsort(index_natsorted(pd_df["CHROM"])))
  renameChrToRoman(pd_df)
  
  fig, ax = plt.subplots()

  renameChrToRoman(df_chr_lengths)
  df_chr_lengths.plot(kind='barh',legend=False, ax=ax, x="chromosome", y="end_position", fontsize=10, figsize=(20,6), edgecolor="black", linewidth=1, color="lightgray", zorder=0, label="Position")

  sns.scatterplot(ax=ax, data=pd_df[pd_df["REF"].isin(["C",'cen'])], x="POS", y="CHROM", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=.9,zorder=1, linewidth=.05, edgecolor="black", legend = False)
  sns.scatterplot(ax=ax, data=pd_df[pd_df["REF"].isin(["G",'cen'])], x="POS", y="CHROM", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=.9,zorder=2, linewidth=.15, edgecolor="black")

  ax.set_xlabel("POSITION", fontsize=10)
  ax.set_ylabel("CHROMOSOME", fontsize=10)
  fig.suptitle(title, fontsize=15)

  legend_texts = plt.legend().get_texts()
  for i in range(len(legend_texts)):
    label_str = legend_texts[i].get_text()
    if label_str in label_dict:
      new_label = label_dict.get(label_str)
      legend_texts[i].set_text(new_label)

  if saveMap:
    plt.savefig('figures/'+title+'.png', transparent=False, dpi=600)
