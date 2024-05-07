import pandas as pd
import numpy as np
import math
import statistics
import logging
import time

logger = logging.getLogger(__name__)

#Loading snakemake rule's variables
input_file = snakemake.input.mztab_input
output = snakemake.output.fasta_output
cutoff = snakemake.params.cutoff

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


df = process_casanovo(input_file)
counter = 0
sequences = df['Casanovo Peptide']
scores = df['Casanovo Score']
with open(output, 'w') as file:
     for sequence, score in zip(sequences, scores):
        header = ">N_" + str(counter) + "_casanovo_" + str(score)
        seq = sequence
        if score >= cutoff:
            counter = counter + 1
            seq = sequence.upper() 
            file.write(header + "\n")
            file.write(seq + "\n")