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


def process_peptideshaker(dbreport_path): 
    """parse peptideshaker report to dataframe (source: https://github.com/denisbeslic/denovopipeline/)

        :param
            dbreport_path: path to the result file of peptideshaker (Extended PSM Report.txt)
        :return
            dbreport_df: dataframe with sequence, related protein and validation
    """
    try:
        with open(dbreport_path) as a:
            dbreport_df = pd.read_csv(dbreport_path, sep="\t")
            validated_list = dbreport_df['Validation'] != 'Not Validated'
            dbreport_df = dbreport_df[validated_list]
            spectrum_title = dbreport_df['Spectrum Title']
            dbreport_df['Index'] = [int(i.split("Index: ")[1].split(",")[0]) for i in spectrum_title]
            db_peptide = dbreport_df['Modified Sequence']
            dbreport_df['Modified Sequence'] = [(i.replace("I", "L").replace("NH2-", "").replace("-COOH", "").replace(
                "C<cmm>", "C").replace("M<ox>", "m").replace("pyro-", "")).replace("N<deam>", "n").replace("Q<deam>",
                                                                                                           "q") for i in
                                                db_peptide]
            db_peptide= dbreport_df['Protein(s)']
            dbreport_df = dbreport_df[['Index', 'Modified Sequence', 'Protein(s)',  'Validation', 'Spectrum Title']]
            dbreport_df = dbreport_df.set_index('Index').sort_index()
            return dbreport_df
    except IOError:
        logger.error(f"Database Search Report from PeptideShaker is not accessible. Make sure it is in {dbreport_path}")
        return pd.DataFrame()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract groundtruth sequences.")
    parser.add_argument("--peptide_shaker_path", type=str, help="Path to Peptideshaker Report.")
    parser.add_argument("--output_file", type=str, help="Output file name.")

    args = parser.parse_args()

    df = process_peptideshaker(args.peptide_shaker_path)

    sequences = df['Modified Sequence']
    proteins = df['Protein(s)']

    with open(args.output_file, "w") as file:
        for sequence, protein in zip(sequences, proteins):
            if "-" in sequence:
                sequence = sequence.replace("-", "")
            header = ">" + protein
            sequence = sequence.replace("<*>", "")
            seq = sequence.upper() 
            file.write(header + "\n")
            file.write(seq + "\n")