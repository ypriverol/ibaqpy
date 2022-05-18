from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

PROTEIN_NAME = 'ProteinName'
GENE_NAME = 'GeneName'
PEPTIDE_SEQUENCE = 'PeptideSequence'
PEPTIDE_CANONICAL = "PeptideCanonical"
PEPTIDE_CHARGE = 'PrecursorCharge'
FRAGMENT_ION = 'FragmentIon'
PRODUCT_CHARGE = 'ProductCharge'
ISOTOPE_LABEL_TYPE = 'IsotopeLabelType'
CONDITION = 'Condition'
BIOREPLICATE = 'BioReplicate'
RUN = 'Run'
FRACTION = 'Fraction'
INTENSITY = 'Intensity'
NORM_INTENSITY = 'NormIntensity'
RT = 'Rt'
REFERENCE = 'Reference'
SAMPLE_ID = 'SampleID'
STUDY_ID = 'StudyID'
SEARCH_ENGINE = 'searchScore'
SCAN = 'Scan'
MBR = 'MatchBetweenRuns'
IBAQ = 'Ibaq'
IBAQ_LOG = 'IbaqLog'
IBAQ_R = "IbaqR"
IBAQ_Q = "IbaqQ"
IBAQ_PPB = 'IbaqPpb'


def remove_contaminants_decoys(dataset: DataFrame, contaminants: str, protein_field=PROTEIN_NAME) -> DataFrame:
    """
    This method reads a file with a list of contaminants and high abudant proteins and
    remove them from the dataset.
    :param dataset: Peptide intensity DataFrame
    :param contaminants_file: contaminants file
    :param protein_field: protein field
    :return: dataset with the filtered proteins
    """
    contaminants_reader = open(contaminants, 'r')
    contaminants_list = contaminants_reader.read().split("\n")
    contaminants_list = [cont for cont in contaminants_list if cont.strip()]

    contaminants_list.append('CONTAMINANT')
    contaminants_list.append('DECOY')
    #cregex = ".*(" + '|'.join(contaminants) + ").*"
    cregex = '|'.join(contaminants_list)
    #for contaminant in contaminants:
        #dataset.drop(index=dataset[dataset[protein_field].str.contains(contaminant)].index, inplace=True)

    return dataset[~dataset[protein_field].str.contains(cregex)]

def plot_distributions(dataset: DataFrame, field: str, class_field: str,
                       file_name: str = None, log2: bool = True, show: bool = False) -> None:
    """
    Print the quantile plot for the dataset
    :param dataset: DataFrame
    :param field: Field that would be use in the dataframe to plot the quantile
    :param class_field: Field to group the quantile into classes
    :param log2: Log the intensity values
    :return:
    """
    pd.set_option('mode.chained_assignment', None)
    normalize = dataset[[field, class_field]]
    if log2:
        normalize[field] = np.log2(normalize[field])
    normalize.dropna(subset=[field], inplace=True)
    data_wide = normalize.pivot(columns=class_field,
                                values=field)
    # plotting multiple density plot
    if show:
        data_wide.plot.kde(figsize=(8, 6), linewidth=2, legend=False)
    if file_name is not None:
        plot = data_wide.plot.kde(figsize=(8, 6), linewidth=2, legend=False)
        plot.figure.savefig(file_name)

    pd.set_option('mode.chained_assignment', 'warn')

def plot_box_plot(dataset: DataFrame, field: str, class_field: str, log2: bool = False, weigth: int = 10,
                  rotation: int = 45, title: str = "", violin: bool = False, file_name: str = None, show: bool = False) -> None:
    """
    Plot a box plot of two values field and classes field
    :param violin: Also add violin on top of box plot
    :param dataset: Dataframe with peptide intensities
    :param field: Intensity field
    :param class_field: class to group the peptides
    :param log2: transform peptide intensities to log scale
    :param weigth: size of the plot
    :param rotation: rotation of the x-axis
    :param title: Title of the box plot
    :return:
    """
    pd.set_option('mode.chained_assignment', None)
    normalized = dataset[[field, class_field]]
    np.seterr(divide='ignore')
    plt.figure(figsize=(weigth, 10))
    if log2:
        normalized[field] = np.log2(normalized[field])

    if violin:
        chart = sns.violinplot(x=class_field, y=field, data=normalized, boxprops=dict(alpha=.3), palette="muted")
    else:
        chart = sns.boxplot(x=class_field, y=field, data=normalized, boxprops=dict(alpha=.3), palette="muted")

    chart.set(title=title)
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation)
    if show:
        plt.show()
    if file_name is not None:
        fig = chart.get_figure()
        fig.savefig(file_name)

    pd.set_option('mode.chained_assignment', 'warn')
