import math
from typing import List, Dict

import click
import pandas as pd
from pandas import DataFrame, Series
from pyopenms import *
import os
import qnorm
import requests

from ibaqpy_commons import remove_contaminants_decoys, PROTEIN_NAME, CONDITION, IBAQ, IBAQ_LOG, IBAQ_PPB, \
    NORM_INTENSITY, SAMPLE_ID, PEPTIDE_CANONICAL, IBAQ_R, IBAQ_Q, plot_distributions, plot_box_plot, GENE_NAME


def uniprot_to_ensembl(uniprot_ids: List[str]) -> Dict:
    """
    Converts a list of uniprot IDs to Ensembl IDs using the Uniprot conversion
    tool.
    Parameters
    ----------
    uniprot_ids : List[str]
        List of IDs to convert.
    Returns
    -------
    DataFrame
        DataFrame with uniprot IDs and corresponding Ensembl IDs.
    """
    chunks = [uniprot_ids[x:x + 100] for x in range(0, len(uniprot_ids), 100)]
    to_df: Dict[str, List] = {}
    deprecated_accessions = []
    for sublist in chunks:
        query = " ".join(sublist)
        url_uniprot = "https://www.uniprot.org/uploadlists/"
        payload = {
            "from": "ACC+ID",
            "to": "GENENAME",
            "format": "tab",
            "query": query,
        }
        r = requests.get(url_uniprot, params=payload)
        print("Request from Uniprot: {}".format(r.url))
        lines = r.text.split("\n")
        for line in lines[1:]:
            try:
                u, e = line.split("\t")
                list_genes = []
                if u in to_df:
                    list_genes = to_df[u]
                list_genes.append(e)
                to_df[u] = list_genes
            except ValueError:
                pass
        # Collect deprecated accessions
        for protein_a in sublist:
            if protein_a not in to_df:
                    deprecated_accessions.append(protein_a)
    # Search deprecated accessions
    if len(deprecated_accessions) > 0:
        deprecated_accessions = list(set(deprecated_accessions))
        for protein_a in deprecated_accessions:
            url_uniprot = "https://www.ebi.ac.uk/proteins/api/uniparc/accession/{}?rfDdtype=Ensembl"
            r = requests.get(url_uniprot.format(protein_a), headers={ "Accept" : "application/json"})
            print("Request from Uniprot: {}".format(url_uniprot.format(protein_a)))
            json_data = r.json()
            gene_list = []
            if "dbReference" in json_data:
                references = json_data["dbReference"]
                for reference in references:
                    if "property" in reference:
                        for property in reference["property"]:
                            if property["type"] == "gene_name":
                                gene_list.append(property["value"])
            gene_list = list(set(gene_list))
            if len(gene_list) > 0:
                to_df[protein_a] = list_genes
    # One protein is maybed associated with multiple genes see example:
    # P62805, in this case, first gene name is selected.
    for protein in to_df:
        sorted_list = to_df[protein]
        sorted_list.sort()
        to_df[protein] = [sorted_list[0]]

    return to_df


def get_canonical_uniprot_accession(accession: str) -> str:
    """
    Protein accession with the structure sp|Q16585|SGCB_HUMAN . To query Uniprot
    the protein ID is used Q16585.
    :param accession: Protein accession with prefix and ID
    :return: return protein id.
    """
    canonical = accession.split("|")[1]
    return canonical.split("-")[0]


def get_gene_accession(accession: str, protein_dic: Dict) -> str:
    accession_list = accession.split(";")
    gene_accessions = []
    for a in accession_list:
        accession = get_canonical_uniprot_accession(a)
        if accession not in protein_dic:
            gene_accessions.append(None)
        else:
            gene_accessions = gene_accessions + protein_dic[accession]
    gene_accessions = list(set(gene_accessions))
    result_id = ""
    remove = False
    for gene in gene_accessions:
        if gene is not None:
            remove = True
    if remove:
        gene_accessions = list(filter(lambda a: a != None, gene_accessions))
    return gene_accessions[0] if len(gene_accessions) == 1 else ";".join(gene_accessions)


def map_protein_gene_expression(res: DataFrame) -> DataFrame:
    """
    Get for a protein dataframe where PROTEIN_NAME is the column of the protein accession
    The method retrive the proteins using the uniprot accession resolver service.
    :param res: protein expression dataframe
    :return: Dictionary with protein accessions with the gene accessions
    """
    protein_list = [e for v in ([x.split(';') for x in res[PROTEIN_NAME].unique()]) for e in v]
    protein_list = [get_canonical_uniprot_accession(x) for x in protein_list]
    print("Number of unique protein accessions: {}".format(len(protein_list)))
    protein_dic = uniprot_to_ensembl(protein_list)

    gene_df = res.copy()
    gene_df[GENE_NAME] = gene_df[PROTEIN_NAME].apply(lambda x: get_gene_accession(x, protein_dic))

    return gene_df


def print_help_msg(command):
    """
    Print the help of the command
    :param command: command
    :return:
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))


def normalize_ibaq(res: DataFrame, quantile: bool) -> DataFrame:
    """
    Normalize the ibaq values using the total ibaq of the sample. The resulted
    ibaq values are then multiplied by 100'000'000 (PRIDE database noramalization)
    for the ibaq ppb and log10 shifted by 10 (ProteomicsDB)
    :param res: Dataframe
    :return:
    """

    res = res[~res.isin([np.nan, np.inf, -np.inf]).any(1)]
    # Get the total intensity for the sample.
    #total_ibaq = res[IBAQ].sum()
    res[IBAQ_R] = res[IBAQ] / res.groupby([CONDITION, SAMPLE_ID])[IBAQ].transform('sum')
    # Normalization method used by Proteomics DB 10 + log10(ibaq/sum(ibaq))
    res[IBAQ_LOG] = res[IBAQ_R].apply(lambda x: 100 + math.log2(x))
    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res[IBAQ_PPB] = res[IBAQ_R].apply(lambda x: x * 100000000)

    if quantile:
        # pivot to have one col per sample
        normalize_df = pd.pivot_table(res, index=[CONDITION, PROTEIN_NAME], columns=SAMPLE_ID, values=IBAQ_PPB,
                                      aggfunc={IBAQ_PPB: np.mean})
        normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(id_vars=[CONDITION, PROTEIN_NAME])
        normalize_df.rename(columns={'value': IBAQ_Q}, inplace=True)
        normalize_df = normalize_df.dropna()
        normalize_df.set_index([PROTEIN_NAME, SAMPLE_ID, CONDITION], inplace=True)
        res.set_index([PROTEIN_NAME, SAMPLE_ID, CONDITION], inplace=True)
        res = pd.concat([res, normalize_df], axis=1)
        res = res.dropna(thresh=1)
        res = res.reset_index()
        res[IBAQ_Q] = res[IBAQ_Q].apply(lambda x: (x) * 100000000)
        print(res.head())

    return res


@click.command()
@click.option("-f", "--fasta", help="Protein database to compute IBAQ values")
@click.option("-p", "--peptides",
              help="Peptide identifications with intensities following the peptide intensity output")
@click.option("-e", "--enzyme", help="Enzyme used during the analysis of the dataset (default: Trypsin)",
              default="Trypsin")
@click.option("-n", "--normalize", help="Normalize IBAQ values using by using the total IBAQ of the experiment",
              is_flag=True)
@click.option("-c", "--contaminants", help="Contaminants protein accession", default="contaminants_ids.tsv")
@click.option("-q", "--quality", help="Applied quality control rules to remove week evidences", is_flag=True)
@click.option("-o", "--output", help="Output file with the proteins and ibaq values")
def ibaq_compute(fasta: str, peptides: str, enzyme: str, normalize: bool, contaminants: str, quality: bool,
                 output: str) -> None:
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

    output_file_prefix = os.path.splitext(output)[0]

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

    # Remove all the proteins with less than 2 peptides per samples
    if quality:
        # Check that peptides are higher than 6 aminoacids
        mask = (data[PEPTIDE_CANONICAL].str.len() > 7)
        dataf = data.loc[mask]
        # Check that proteins contains more than 2 unique peptides per proteins.
        dataf = dataf.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].count()
        dataf = dataf[dataf > 1]
        data.set_index([PROTEIN_NAME, SAMPLE_ID, CONDITION], inplace=True)
        data = data.loc[dataf.index, :]
        data = data.reset_index()
        print(data)

    res = pd.DataFrame(data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()).apply(
        get_average_nr_peptides_unique_bygroup, 1)
    res = res.sort_values(ascending=False)
    res = res.reset_index()
    res = res.rename(columns={0: IBAQ})

    # If samples are higher than 3, the protein MUST be identified in atl east 2 samples.
    num_samples = len(res[SAMPLE_ID].unique())
    if quality and num_samples > 3:
        dataf = res.groupby([PROTEIN_NAME, CONDITION])[SAMPLE_ID].count()
        dataf = dataf[dataf > 1]
        dataf = dataf.reset_index()
        print("Number of proteins: " + str(len(res[PROTEIN_NAME].unique())))
        res = res[res[PROTEIN_NAME].isin(dataf[PROTEIN_NAME].values)]
        print("Number of proteins in more than 2 Samples: " + str(len(res[PROTEIN_NAME].unique())))
        print(res)

    res = remove_contaminants_decoys(res, contaminants, protein_field=PROTEIN_NAME)
    if normalize:
        res = normalize_ibaq(res, quantile=True)

    gene_expression_df = map_protein_gene_expression(res)
    gene_expression_df = gene_expression_df.groupby([GENE_NAME, SAMPLE_ID, CONDITION], as_index=True)[[IBAQ, IBAQ_LOG, IBAQ_PPB]].mean()
    gene_expression_df = gene_expression_df.reset_index()

    num_samples = len(gene_expression_df[SAMPLE_ID].unique())
    if quality and num_samples > 3:
        dataf = gene_expression_df.groupby([GENE_NAME, CONDITION])[SAMPLE_ID].count()
        dataf = dataf[dataf > 1]
        dataf = dataf.reset_index()
        print("Number of genes: " + str(len(gene_expression_df[GENE_NAME].unique())))
        gene_expression_df = gene_expression_df[gene_expression_df[GENE_NAME].isin(dataf[GENE_NAME].values)]
        print("Number of genes in more than 2 Samples: " + str(len(gene_expression_df[GENE_NAME].unique())))
        print(gene_expression_df)

    print("Printing the quantile IBAQ distributions")
    plot_distributions(gene_expression_df, IBAQ_PPB, SAMPLE_ID, log2=False,
                       file_name=output_file_prefix + "-{}-{}.pdf".format("gene-distribution", "IBAQ_PPB"), show=False)
    plot_box_plot(gene_expression_df, IBAQ_PPB, SAMPLE_ID, log2=True, title="Gene IBAQ Quantile Distribution", violin=False,
                  file_name=output_file_prefix + "-{}-{}.pdf".format("gene-boxplot", "IBAQ_PPB"), show=False)

    print("Printing the rIBAQ distributions")
    plot_distributions(res, IBAQ_PPB, SAMPLE_ID, log2=True,
                       file_name=output_file_prefix + "-{}-{}.pdf".format("distribution", "IBAQ_PPB"), show=False)
    plot_box_plot(res, IBAQ_PPB, SAMPLE_ID, log2=True, title="Original IbaqLOG (no normalization)", violin=False,
                  file_name=output_file_prefix + "-{}-{}.pdf".format("boxplot", "IBAQ_PPB"), show=False)

    res.to_csv(output, index=False)
    gene_expression_df.to_csv(output_file_prefix + "-gene.tsv", index=False)


if __name__ == '__main__':
    ibaq_compute()
