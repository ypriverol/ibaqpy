import math

import click
import pandas as pd
from pandas import DataFrame, Series
from pyopenms import *

from ibaqpy_commons import remove_contaminants_decoys, PROTEIN_NAME, CONDITION, IBAQ, IBAQ_LOG, IBAQ_PPB, \
    NORM_INTENSITY, SAMPLE_ID, PEPTIDE_CANONICAL


def print_help_msg(command):
    """
    Print the help of the command
    :param command: command
    :return:
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))

def normalize_ibaq(res: DataFrame) -> DataFrame:
    """
    Normalize the ibaq values using the total ibaq of the sample. The resulted
    ibaq values are then multiplied by 100'000'000 (PRIDE database noramalization)
    for the ibaq ppb and log10 shifted by 10 (ProteomicsDB)
    :param res: Dataframe
    :return:
    """

    res = res[~res.isin([np.nan, np.inf, -np.inf]).any(1)]
    # Get the total intensity for the sample.
    total_ibaq = res[IBAQ].sum()
    # Normalization method used by Proteomics DB 10 + log10(ibaq/sum(ibaq))
    res[IBAQ_LOG] = res[IBAQ].apply(lambda x: 100 + math.log2(x / total_ibaq))
    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res[IBAQ_PPB] = res[IBAQ].apply(lambda x: (x / total_ibaq) * 100000000)
    return res

@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option("-p", "--peptides", help="Peptide identifications with intensities following the peptide intensity output")
@click.option("-e", "--enzyme", help="Enzyme used during the analysis of the dataset (default: Trypsin)",
              default="Trypsin")
@click.option("-n", "--normalize", help="Normalize IBAQ values using by using the total IBAQ of the experiment",
              is_flag=True)
@click.option("-c", "--contaminants", help="Contaminants protein accession", default="contaminants_ids.tsv")
@click.option("-q", "--quality", help="Applied quality control rules to remove week evidences", is_flag = True)
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
def ibaq_compute(fasta: str, peptides: str, enzyme: str, normalize: bool, contaminants: str, quality: bool, output: str) -> None:
    """
    This command computes the IBAQ values for a file output of peptides with the format described in
    peptide_contaminants_file_generation.py.

    quality control is based on two main rules:
    - Proteins MUST contain at least two unique peptides to be considered identified
    - Identified peptides MUST contains at least 8 AAs.

    :param fasta: Fasta file used to perform the peptide identification
    :param peptides: Peptide intensity file
    :param enzyme: Enzyme used to digest the protein sample
    :param normalize: use some basic normalization steps.
    :param contaminants_file: Contaminant data file
    :param output: output format containing the ibaq values.
    :return:
    """
    if peptides is None or fasta is None:
        print_help_msg(ibaq_compute)
        exit(1)

    fasta_proteins = list()  # type: list[FASTAEntry]
    FASTAFile().load(fasta, fasta_proteins)
    uniquepepcounts = dict()  # type: dict[str, int]
    MINLEN = 6
    MAXLEN = 30
    ENZYMENAME = enzyme
    digestor = ProteaseDigestion()
    digestor.setEnzyme(ENZYMENAME)

    def get_average_nr_peptides_unique_bygroup(pdrow: Series) -> Series:
        """
        Get the average intensity for protein gorups
        :param pdrow: peptide row
        :return: average intensity
        """
        proteins = pdrow.name[0].split(';')
        summ = 0
        for prot in proteins:
            summ += uniquepepcounts[prot]
        return pdrow.NormIntensity / (summ / len(proteins))

    for entry in fasta_proteins:
        digest = list()  # type: list[str]
        digestor.digest(AASequence().fromString(entry.sequence), digest, MINLEN, MAXLEN)
        digestuniq = set(digest)
        uniquepepcounts[entry.identifier] = len(digestuniq)

    data = pd.read_csv(peptides, sep="\t")
    print(data.head())
    # next line assumes unique peptides only (at least per indistinguishable group)

    dataf = None
    if quality:
      # Check that peptides are higher than 6 aminoacids
      mask = (data[PEPTIDE_CANONICAL].str.len() > 7)
      dataf = data.loc[mask]
      # Check that proteins contains more than 2 unique peptides per proteins.
      dataf = dataf.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].count()
      dataf = dataf[dataf > 1]

    res = pd.DataFrame(data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()).apply(get_average_nr_peptides_unique_bygroup, 1)
    res = res.sort_values(ascending=False)
    res = res.reset_index()
    res = res.rename(columns={0: IBAQ})

    res = remove_contaminants_decoys(res, contaminants, protein_field=PROTEIN_NAME)
    if normalize:
        res = normalize_ibaq(res)

    if quality and dataf is not None:
        res.set_index([PROTEIN_NAME, SAMPLE_ID, CONDITION], inplace=True)
        res = res.loc[dataf.index,:]
        res = res.reset_index()
        print(res)

    res.to_csv(output, index=False)

if __name__ == '__main__':
    ibaq_compute()
