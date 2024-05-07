import pandas as pd
import numpy as np
import math
import statistics
import logging
import time
import os
import argparse

logger = logging.getLogger(__name__)


def process_casanovo(casanovo_path):
    """
    ###### From https://github.com/denisbeslic/denovopipeline and modified by AL ######
    parse de novo results from Casanovo to dataframe
            :param
                casanovo_path: path to the mztab file of Casanovo
            :return
                pnovo_df: dataframe with Score, Peptide, AAScore of Casanovo
    """
    try:
        with open(casanovo_path) as f:

            # Identify header
            for i, line in enumerate(f):
                if line.startswith("PSH"):
                    casanovo_header = i
                    break

            casanovo_df = pd.read_csv(casanovo_path, sep="	", header=casanovo_header)
            casanovo_peptide = casanovo_df['sequence'].tolist()
            casanovo_df['Casanovo Peptide'] = [str(i).replace('+15.995', '').replace('+0.984', '').replace(' ',
                                        '').replace('+57.021', '').replace('+42.011', '').replace('-17.027',
                                        '').replace('+43.006', '') for i in casanovo_peptide]
            casanovo_df['Casanovo Score'] = casanovo_df['search_engine_score[1]']
            casanovo_df['Casanovo Score'] = [i * 100 for i in casanovo_df['Casanovo Score'].tolist()]
            casanovo_aascore = casanovo_df['opt_ms_run[1]_aa_scores'].tolist()
            casanovo_aascore= [i.replace(',', ' ') for i in casanovo_aascore]
            casanovo_aascore = [i.split() for i in casanovo_aascore]
            casanovo_aascore = [[str(float(j) * 100) for j in i] for i in casanovo_aascore]
            casanovo_df['Casanovo aaScore'] = [" ".join(i) for i in casanovo_aascore]
            casanovo_df = casanovo_df[['Casanovo Peptide', 'Casanovo Score', 'Casanovo aaScore']]
            return casanovo_df
    except IOError:
        logger.warning(f"Casanovo results not accessible. Make sure they are placed in {casanovo_path}")
        return pd.DataFrame()

def process_smsnet(smsnet_path):
    """
    ###### From https://github.com/denisbeslic/denovopipeline and modified by AL ######
    parse de novo results from SMSNet to dataframe

            :param
                smsnet_path: path to the result file of Novor
            :return
                smsnet_df: dataframe with Score, Peptide, AAScore of SMSNet
        """
    try:
        with open(smsnet_path) as f, open(
                smsnet_path + '_prob') as g:  # change _rescore and _prob to switch between rescoring and real
            smsnet_peptide = pd.Series([line.rstrip() for line in f])
            peptide_list = [x.replace(" ", "").replace("I", "L").replace("<s>","").replace("<unk>","") for x in
                            smsnet_peptide]
            smsnet_peptide = pd.DataFrame(peptide_list)
            aa_score = g.readlines()
            aa_score = [i.strip().split(' ') for i in aa_score]
            score_sum = []
            for i in range(len(aa_score)):
                if not aa_score[i] == ['']:
                    for j in range(len(aa_score[i])):
                        aa_score[i][j] = float(np.exp(float(aa_score[i][j])) * 100)
                else:
                    aa_score[i] = [0]
             # Compute peptide score and set as SMSNet score, the mean of aaScores
            for i in range(len(aa_score)):
                if not aa_score[i] == [0]:
                    score_sum.append(statistics.mean(aa_score[i]))
                else:
                    score_sum.append(0)
            df = pd.DataFrame({'aaScore': aa_score, 'Peptide Score': score_sum})
            smsnet_df = pd.concat([smsnet_peptide, df], axis=1)
            smsnet_df.columns = ['SMSNet Peptide', 'SMSNet aaScore', 'SMSNet Score']
            clist = ['SMSNet Score', 'SMSNet Peptide', 'SMSNet aaScore']
            smsnet_df = smsnet_df[clist]
            smsnet_df.index = range(1, len(smsnet_df) + 1)

            smsnet_aascore = smsnet_df['SMSNet aaScore'].tolist()

            # remove peptide predictions which look like "<s><s><ink>SSSSLASSS"
            smsnet_df['SMSNet aaScore'] = [str(i).replace(' ', '').replace(',', ' '
                                            ).replace('[', '').replace(']', '')
                                           for i in
                                           smsnet_aascore]
            smsnet_df = smsnet_df[['SMSNet Peptide', 'SMSNet Score', 'SMSNet aaScore']]

            return smsnet_df
    except IOError:
        logger.warning(f"SMSNet results not accessible. Make sure they are placed in {smsnet_path}")
        return pd.DataFrame()

def naive_search(df_1_, df_2_, pNovo = False, merged = False):
  """
  Naive search between peptide sequences

  Parameters:
      df_1 (data frame): first data frame contains peptide sequences and scores
      df_2 (data frame): second data frame contains peptide sequences and scores
      pNovo (bool): check if the input data frame is from pNovo results as pNovo has only Score
      merged (bool): check if the first dataframe is already concatenation of others!

  Returns:
    output_def (data frame): results summary (only shared sequences)
    matches (int): total number of shared peptides
  """

  df_1 = df_1_
  df_2 = df_2_
  summary_1 = df_1_.columns[0].split()
  summary_2 = df_2_.columns[0].split()
  complete_matches = 0

  df_1_alg = ""
  df_2_alg = "

  #Algorithms names
  df_1_alg = df_1.columns[0].split()
  df_1_alg_split = df_1_alg[0]

  df_2_alg = df_2.columns[0].split()
  df_2_alg_split = df_2_alg[0]

  #Store tool headers
  headers_df_1 = list()
  headers_df_2 = list()

  headers_df_1.append(df_1.columns[0])
  headers_df_2.append(df_2.columns[0])

  #Change headers names for simple search
  df_1.rename(columns={headers_df_1[0]: 'Peptide'}, inplace=True)
  df_2.rename(columns={headers_df_2[0]: 'Peptide'}, inplace=True)

  #Improve performance of counting
  peptides_set_1 = set(df_1['Peptide'].str.upper().str.strip())
  peptides_set_2 = set(df_2['Peptide'].str.upper().str.strip())

  completely_shared_peptides = peptides_set_1.intersection(peptides_set_2)

  #Count unique elements to each set
  unique_peptides_to_set1 = peptides_set_1.difference(peptides_set_2)
  unique_peptides_to_set2 = peptides_set_2.difference(peptides_set_1)

  #Find elements unique to either set
  unique_to_either_algorithms = unique_peptides_to_set1.symmetric_difference(unique_peptides_to_set2) #Union of unique elements (not shared between two sets)

  #Compute union for Jaccard Index
  union = unique_peptides_to_set1.union(unique_peptides_to_set2)

  Jaccard_index = (len(completely_shared_peptides) / len(union))
  Jaccard_distance = 1 - Jaccard_index

  #Create new dataframe with the old headers names
  peptide_name = 'Peptide_'+ df_1_alg_split + '_' + df_2_alg_split
  output_df = pd.DataFrame(columns=['Peptide'])

  for shared_peptide in completely_shared_peptides:
      output_df.loc[len(output_df)] = [shared_peptide]

  matches_ = len(completely_shared_peptides)
  print("Number of complete matches: ", len(completely_shared_peptides))
  print("Number of unique peptides to both sets: ", len(unique_to_either_algorithms))
  print("Similarity ration (Jaccard Index): ", Jaccard_index*100, "%")
  print("Dissimilarity ratio (Jaccard Distance): ", Jaccard_distance*100, "%")
  print("Number of unique peptides to set 1: ", len(unique_peptides_to_set1), "/", len(df_1_))
  print("Number of unique peptides to set 2: ", len(unique_peptides_to_set2), "/", len(df_2_))

  return output_df, matches_

def unshared_peptides(df_1_, df_2_, pNovo = False, merged = False):
  """
  Naive search between peptide sequences

  Parameters:
      df_1 (data frame): first data frame contains peptide sequences and scores
      df_2 (data frame): second data frame contains peptide sequences and scores
      pNovo (bool): check if the input data frame is from pNovo results as pNovo has only Score
      merged (bool): check if the first dataframe is already concatenation of others!

  Returns:
    output_def (data frame): results summary (only not shared sequences)
    matches (int): total number of not shared peptides
  """

  df_1 = df_1_
  df_2 = df_2_
  summary_1 = df_1_.columns[0].split()
  summary_2 = df_2_.columns[0].split()
  complete_matches = 0


  df_1_alg = ""
  df_2_alg = ""

  df_1_alg = df_1.columns[0].split()
  df_1_alg_split = df_1_alg[0]

  df_2_alg = df_2.columns[0].split()
  df_2_alg_split = df_2_alg[0]


  #Store tool headers
  headers_df_1 = list()
  headers_df_2 = list()

  headers_df_1.append(df_1.columns[0])
  headers_df_2.append(df_2.columns[0])

  #Change headers names for simple search
  df_1.rename(columns={headers_df_1[0]: 'Peptide'}, inplace=True)
  df_2.rename(columns={headers_df_2[0]: 'Peptide'}, inplace=True)

  #Improve performance of counting
  peptides_set_1 = set(df_1['Peptide'].str.upper().str.strip())
  peptides_set_2 = set(df_2['Peptide'].str.upper().str.strip())

  completely_shared_peptides = peptides_set_1.intersection(peptides_set_2)

  #Count unique elements to each set
  unique_peptides_to_set1 = peptides_set_1.difference(peptides_set_2)
  unique_peptides_to_set2 = peptides_set_2.difference(peptides_set_1)

  #Find elements unique to either set
  unique_to_either_algorithms = unique_peptides_to_set1.symmetric_difference(unique_peptides_to_set2) #Union of unique elements (not shared between two sets)

  #Find unshared elements of set_1
  unique_to_set_1 = unique_peptides_to_set1.difference(unique_peptides_to_set2)
  #Compute union for Jaccard Index
  union = unique_peptides_to_set1.union(unique_peptides_to_set2)

  #Create new dataframe with the old headers names
  peptide_name = 'Peptide_'+ df_1_alg_split + '_' + df_2_alg_split
  output_df = pd.DataFrame(columns=['Peptide'])

  for not_shared_peptide in unique_to_set_1:
      output_df.loc[len(output_df)] = [not_shared_peptide]

  non_matches_ = len(unique_to_set_1)

  return output_df, non_matches_

def generate_fasta(sample_name, casanovo_path, smsnet_path, output_path, score_cutoff):

    """
    Generate fasta files of the merged process of casanovo and SMSNet

        :param
            sample_name: name of the sample, will be added to the headers and to the output fasta file
            casanovo_path: path to the result of casanovo peptide sequences
            smsnet_path: path to the result of smsnet peptide sequences
            output_path: path to the output directory
            score_cutoff: score cutoff of the de-novo peptide sequences
        :output
            fasta_file: the resulted fasta file with the header:
            > (sample_name) + (shared/not_shared) + casanovo_score and/or smsnet_score
            AAALLLYYNNMAMMAA....
         :return
            None
    """

    print(f'[Info] Reading raw peptides predictions......\n')
    print(f' Casanovo: {casanovo_path}')
    print(f' SMSNet: {smsnet_path}')

    start_time = time.time()

    df_casanovo = process_casanovo(casanovo_path)
    df_smsnet = process_smsnet(smsnet_path)

    shared_sequences = dict()
    non_shared_casanovo = dict()
    non_shared_smsnet = dict()
    fasta_dictionary = dict()

    print(f'[Info] Finsihed loading sequnces! \n')

    print(f'[Info] Jaccard Index computations started \n')
    print(f'[Info] Esitmating number of shared sequences.... \n')
    shared_casanovo_smsnet, casanovo_smsnet = naive_search(df_casanovo, df_smsnet, False, False)

    print(f'[Info] Detected {casanovo_smsnet} shared sequences. \n')

    print(f'[Info] Converting all shared sequences to dicionary .... \n')

    for shared_peptide in shared_casanovo_smsnet['Peptide']:
      casanovo_score = df_casanovo[df_casanovo['Peptide'].str.upper().str.strip() == shared_peptide][df_casanovo.columns[1]].iloc[0]
      if casanovo_score < 0 : casanovo_score *= -1
      smsnet_score = df_smsnet[df_smsnet['Peptide'].str.upper().str.strip() == shared_peptide][df_smsnet.columns[1]].iloc[0]
      if smsnet_score < 0 : smsnet_score *= -1
      if casanovo_score >= score_cutoff or smsnet_score >= score_cutoff:
        sequence_header = f'>{sample_name}|shared|[casanovoScore: {casanovo_score}]|[smsnetScore: {smsnet_score}]'
        shared_sequences[sequence_header] = shared_peptide
      else:
        continue

    print(f'[Info] The dicionary of shared sequences is created with {len(shared_sequences)} observed peptides. \n')


    print(f'[Info] Converting all casanovo non-shared sequences to dicionary .... \n')

    nonshared_casanovo_smsnet, casanovo_smsnet_non_matches = unshared_peptides(df_casanovo, df_smsnet, False, False)
    merged_df_casanovo = pd.merge(nonshared_casanovo_smsnet, df_casanovo, on='Peptide', how='left')
    duplicate_peptides = merged_df_casanovo[merged_df_casanovo.duplicated(subset='Peptide', keep=False)]

    # Sort the DataFrame by 'Casanovo Score' in descending order
    merged_df_casanovo_sorted = merged_df_casanovo.sort_values(by='Casanovo Score', ascending=False)
    duplicate_peptides = merged_df_casanovo_sorted[merged_df_casanovo_sorted.duplicated(subset='Peptide', keep=False)]

    # Drop duplicates based on 'Peptide', keeping the one with the highest 'Casanovo Score'
    filtered_df_casanovo = merged_df_casanovo_sorted.drop_duplicates(subset='Peptide', keep='first')
    filtered_df_casanovo = filtered_df_casanovo.reset_index(drop=True)

    skipped_casanovo = 0
    for i in range(len(filtered_df_casanovo)):
      peptide = filtered_df_casanovo['Peptide'][i]
      casanovo_score = filtered_df_casanovo['Casanovo Score'][i]
      if casanovo_score < 0 : casanovo_score *= -1
      if casanovo_score >= score_cutoff:
        sequence_header = f'>{sample_name}|casanovo-non-shared|[casanovoScore: {casanovo_score}]'
        non_shared_casanovo[sequence_header] = peptide
      else:
        skipped_casanovo += 1
        continue

    print(f'[Info] The dicionary of casanovo non-shared sequences is created with {len(non_shared_casanovo)} observed peptides. \n')

    print(f'[Info] Converting all smsnet non-shared sequences to dicionary .... \n')

    nonshared_smsnet_casanovo, casanovo_smsnet_non_matches = unshared_peptides(df_smsnet, df_casanovo, False, False)
    merged_df_smsnet = pd.merge(nonshared_smsnet_casanovo, df_smsnet, on='Peptide', how='left')
    duplicate_peptides = merged_df_smsnet[merged_df_smsnet.duplicated(subset='Peptide', keep=False)]

    # Sort the DataFrame by 'SMSNet Score' in descending order
    merged_df_smsnet_sorted = merged_df_smsnet.sort_values(by='SMSNet Score', ascending=False)
    duplicate_peptides = merged_df_smsnet_sorted[merged_df_smsnet_sorted.duplicated(subset='Peptide', keep=False)]

    # Drop duplicates based on 'Peptide', keeping the one with the highest 'SMSNet Score'
    filtered_df_smsnet = merged_df_smsnet_sorted.drop_duplicates(subset='Peptide', keep='first')
    filtered_df_smsnet = filtered_df_smsnet.reset_index(drop=True)

    skipped_smsnet = 0
    for i in range(len(filtered_df_smsnet)):
      peptide = filtered_df_smsnet['Peptide'][i]
      smsnet_score = filtered_df_smsnet['SMSNet Score'][i]
      if smsnet_score >= score_cutoff:
        sequence_header = f'>{sample_name}|smsnet-non-shared|[smsnetScore: {smsnet_score}]'
        non_shared_smsnet[sequence_header] = peptide
      else:
        skipped_smsnet += 1
        continue

    print(f'[Info] The dicionary of smsnet non-shared sequences is created with {len(non_shared_smsnet)} observed peptides. \n')

    print(f'[Info] Combining all three dicionaries into one .... \n')
    fasta_dictionary.update(shared_sequences)
    fasta_dictionary.update(non_shared_casanovo)
    fasta_dictionary.update(non_shared_smsnet)
    print(f'[Info] Combination is finished with {len(fasta_dictionary)} observed peptides. \n')

    #Removing all empty values, this is because an writing error of SMSNet
    empty_nan = [key for key, value in fasta_dictionary.items() if value == '']
    for key in empty_nan:
        fasta_dictionary.pop(key)
    
    print("[Info] Total number of casanovo skipped peptides: ", skipped_casanovo)
    print("[Info] Total number of smsnet skipped peptides: ", skipped_smsnet)
    print(f'[Info] Writing sequences to the fasta file {output_path}/{sample_name}.fasta.... \n')
    #output = f"{output_path}/{sample_name}.fasta"
    output = output_path

    with open(output, 'w') as output_file:
        for header, sequence in fasta_dictionary.items():
            output_file.write(header + '\n')  # Write the header line
            output_file.write(sequence + '\n')  # Write the sequence
    
    print("[Info] Total number of generated peptides: ", len(fasta_dictionary))
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Elapsed time for the sample {sample_name}: {elapsed_time:.4f} seconds")

    return None

#casanovo_path = snakemake.input.casanovo_path
#smsnet_dir = snakemake.input.smsnet_dir_path
#output_path = snakemake.output.output_dir
#cutoff = snakemake.params.cutoff

#sample_name = smsnet_dir[0].split("/")[-1]
#smsnet_path = smsnet_dir[0] + "/" + sample_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine Casanovo and SMSNet Peptides with Set based Algorithms.")
    parser.add_argument("--sample_name", type=str, help="Name of the sample")
    parser.add_argument("--casanovo_path", type=str, help="Path to Casanovo")
    parser.add_argument("--smsnet_path", type=str, help="Path to SMSNet directory")
    parser.add_argument("--output_path", type=str, help="Path to output directory")
    parser.add_argument("--cutoff", type=float, help="Cutoff value")

    args = parser.parse_args()
    
    generate_fasta(args.sample_name, args.casanovo_path, args.smsnet_path, args.output_path, args.cutoff)

###################################### Metaproteomics 
# samples = [
#     "L00838_S05", "L00842_S02", "L00846_S06", "L00850_S23", 
#     "L00854_S16", "L00858_S25", "L00862_S07", "L00866_S09", 
#     "L00870_S04", "L00874_S30", "L00880_S32", "L00884_S24", 
#     "L00840_S11", "L00844_S08", "L00848_S22", "L00852_S14", 
#     "L00856_S21", "L00860_S03", "L00864_S15", "L00868_S26", 
#     "L00872_S01", "L00876_S31", "L00882_S33", "L00887_S29"
# ]

# for sample_name in samples:
#     casanovo_path = "/home/lutfia/scratch/casanovo/" + sample_name + ".mztab"
#     smsnet_path = "/home/lutfia/scratch/SMSNet/output/" + sample_name
#     output_path = "/home/lutfia/scratch/megaX/build/main/Meta_proteomics_Data/denovo_peptides_casanovo_smsnet/"
#     generate_fasta(sample_name, casanovo_path, smsnet_path, output_path, 60)
