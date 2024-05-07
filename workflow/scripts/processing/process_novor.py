from re import I
from turtle import color
import numpy as np
import pandas as pd
import statistics
import math
import logging
import os 
import argparse
logger = logging.getLogger(__name__)


def process_novor(novor_path):
    """parse de nov results from novor to dataframe (source: https://github.com/denisbeslic/denovopipeline/)

        :param
            novor_path: path to the result file of Novor
        :return
            novor_df: dataframe with ScanNum, Score, Peptide, AAScore of Novor
    """
    try:
        with open(novor_path) as f:

            # jump over 17 rows which are the header of a Novor result file
            novor_df = pd.read_csv(novor_path, sep=",", header=17)

            # change single_score from "14-52-12.." to "14 52 12.."
            single_score = novor_df[' aaScore'].tolist()
            for i in range(len(single_score)):
                single_score[i] = single_score[i].replace(" ", "").replace("-", " ")
            novor_df[' aaScore'] = single_score
            novor_df = novor_df[[' scanNum', ' score', ' peptide', ' aaScore']]
            novor_df = novor_df.set_index(' scanNum')
            novor_df.columns = ['Novor Score', 'Novor Peptide', 'Novor aaScore']
            # Replace the specific annotation of Novor for Modifications
            novor_peptide = novor_df['Novor Peptide'].tolist()
            novor_df['Novor Peptide'] = [i.replace('M(0)', 'm').replace('Q(2)', 'q').replace('N(1)', 'n').replace(' ',
                                        '').replace('C(3)', 'C').replace('Q(0)','q').replace('M(2)','m').replace('K(1)','K').replace('C(2)', 'C').replace('R(0)', 'R') for i in novor_peptide]
            return novor_df
    except IOError:
        logger.warning(f"Novor results not accessible. Make sure they are placed in {novor_path}")
        return pd.DataFrame()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Novor sequences.")
    parser.add_argument("--novor_path", type=str, help="Path to Novor file.")
    parser.add_argument("--output_file", type=str, help="Output file name.")
    parser.add_argument("--cutoff", type=float, help="Cutoff value")

    args = parser.parse_args()
    
    df = process_novor(args.novor_path)
    sequences = df['Novor Peptide']
    scores = df['Novor Score']

    with open(args.output_file, "w") as file:
        count = 0
        for sequence, score in zip(sequences, scores):
            count = count + 1
            if "-" in sequence:
                sequence = sequence.replace("-", "")
            if "(" in sequence:
                sequence = sequence.replace("(", "")
            if ")" in sequence:
                sequence = sequence.replace(")", "")
            header = ">Protein_" + str(score) + '_' + str(count)
            if score >= args.cutoff:
                seq = sequence.upper() 
                file.write(header + "\n")
                file.write(seq + "\n")