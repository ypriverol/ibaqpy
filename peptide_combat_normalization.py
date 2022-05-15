from typing import List

from combat.pycombat import pycombat
import pandas as pd
import matplotlib.pyplot as plt
import click
import glob
import numpy as np
import seaborn as sns

# # prepare data
# # the datasets are dataframes where:
#     # the indexes correspond to the gene names
#     # the column names correspond to the sample names
# # Any number (>=2) of datasets can be treated
# dataset_1 = pd.read_pickle("data/GSE18520.pickle") # datasets can also be stored in csv, tsv, etc files
# dataset_2 = pd.read_pickle("data/GSE66957.pickle")
# dataset_3 = pd.read_pickle("data/GSE69428.pickle")
#
# # we merge all the datasets into one, by keeping the common genes only
# df_expression = pd.concat([dataset_1,dataset_2,dataset_3],join="inner",axis=1)
#
# # plot raw data
# plt.boxplot(df_expression.transpose())
# plt.show()
#
# batch = []
# datasets = [dataset_1,dataset_2,dataset_3]
# for j in range(len(datasets)):
#     batch.extend([j for _ in range(len(datasets[j].columns))])
#
# # run pyComBat
# df_corrected = pycombat(df_expression,batch)
#
# # visualise results
# plt.boxplot(df_corrected.transpose())
# plt.show()
from pandas import DataFrame

from ibaqpy_commons import PEPTIDE_SEQUENCE, CONDITION, INTENSITY, SAMPLE_ID
import os

def print_help_msg(command) -> None:
    """
    Print help information
    :param command: command to print helps
    :return: print help
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))




@click.command()
@click.option("--input", help="Folder with all the peptide tables")
@click.option("--pattern", help="Pattern extension of the file .tsv, *Intensities-Norm.tsv")
@click.option("--verbose",
              help="Print addition information about the distributions of the intensities, number of peptides remove after normalization, etc.",
              is_flag=True)
def peptide_combat_normalization(input: str, pattern: str, verbose: str) -> None:

    if pattern is None:
        print_help_msg(peptide_combat_normalization)
        exit(1)
    if input is None:
        input = "."
    file_list = [f for f in os.listdir(input) if pattern in f]
    dataframes = list()
    for file in file_list:
        df = pd.read_csv(input + "/" + file, sep="\t")
        dataframes.append(df)

    batch = []
    for j in range(len(file_list)):
       batch.extend([j for _ in range(len(dataframes[j].columns))])

    datasets = pd.concat(dataframes,join="inner",axis=1)
    datasets = datasets.fillna(25)
    df_corrected = pycombat(datasets, batch)
    plt.figure(figsize=(15, 10))
    melted = df_corrected.melt()
    melted["value"] = np.log2(melted["value"])
    chart = sns.boxplot(x=SAMPLE_ID, y="value", data=melted, boxprops=dict(alpha=.3), palette="muted")

    chart.set_xticklabels(chart.get_xticklabels(), rotation=45)
    plt.show()

    print(file_list)



if __name__ == '__main__':
    peptide_combat_normalization()