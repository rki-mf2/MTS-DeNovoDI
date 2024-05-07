import pandas as pd
import numpy as np
import math
import statistics
import logging
import time
import os

logger = logging.getLogger(__name__)


#Loading snakemake rule's variables
input_file = snakemake.input.input_file
output = snakemake.output.fasta_output
cutoff = snakemake.params.cutoff

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



input_file = snakemake.input.input_file
output = snakemake.output.fasta_output
cutoff = snakemake.params.cutoff

sample_name = input_file[0].split("/")[-1]
smsnet_path = input_file[0] + "/" + sample_name
df = process_smsnet(smsnet_path)
counter = 0
sequences = df['SMSNet Peptide']
scores = df['SMSNet Score']
with open(output, 'w') as file:
     for sequence, score in zip(sequences, scores):
        header = ">N_" + str(counter) + "_smsnet_" + str(score)
        seq = sequence
        if score >= cutoff:
            counter = counter + 1
            seq = sequence.upper() 
            file.write(header + "\n")
            file.write(seq + "\n")