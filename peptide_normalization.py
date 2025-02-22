import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas import DataFrame
import qnorm
import re
import os
import secrets

from ibaqpy_commons import remove_contaminants_decoys, INTENSITY, SAMPLE_ID, NORM_INTENSITY, \
    PEPTIDE_SEQUENCE, CONDITION, PEPTIDE_CHARGE, FRACTION, RUN, BIOREPLICATE, RT, PEPTIDE_CANONICAL, SEARCH_ENGINE, \
    PROTEIN_NAME, plot_distributions, plot_box_plot


def print_dataset_size(dataset: DataFrame, message: str, verbose: bool) -> None:
    if verbose:
        print(message + str(len(dataset.index)))


def print_help_msg(command) -> None:
    """
    Print help information
    :param command: command to print helps
    :return: print help
    """
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))

def remove_outliers_iqr(dataset: DataFrame):
    """
    This method removes outliers from the dataframe inplace, the variable used for the outlier removal is Intensity
    :param dataset: Peptide dataframe
    :return: None
    """
    Q1 = dataset[INTENSITY].quantile(0.25)
    Q3 = dataset[INTENSITY].quantile(0.75)
    IQR = Q3 - Q1

    dataset.query('(@Q1 - 1.5 * @IQR) <= Intensity <= (@Q3 + 1.5 * @IQR)', inplace=True)

def remove_missing_values(normalize_df: DataFrame, ratio: float = 0.3) -> DataFrame:
    """
    Remove missing values if the peptide do not have values in the most of the samples
    :param normalize_df: data frame with the data
    :param ratio: ratio of samples without intensity values.
    :return:
    """
    n_samples = len(normalize_df.columns)
    normalize_df = normalize_df.dropna(thresh=round(n_samples * ratio))
    return normalize_df

def get_canonical_peptide(peptide_sequence: str)-> str:
    """
    This function returns a peptide sequence without the modification information
    :param peptide_sequence: peptide sequence with mods
    :return: peptide sequence
    """
    clean_peptide = re.sub("[\(\[].*?[\)\]]", "", peptide_sequence)
    clean_peptide = clean_peptide.replace(".","")
    return clean_peptide

def intensity_normalization(dataset: DataFrame, field: str, class_field: str = "all", scaling_method: str = "msstats", indisk: bool = False) -> DataFrame:

    ## TODO add imputation and/or removal to those two norm strategies
    if scaling_method == 'msstats':
        g = dataset.groupby(['Run', 'Fraction'])[INTENSITY].apply(np.median)
        g.name = 'RunMedian'
        dataset = dataset.join(g, on=['Run', 'Fraction'])
        # TODO might be quicker with transform but I could not make it work in short time for multidim. groupby
        #dataset['RunMedian'] = dataset[INTENSITY].groupby(dataset[['Run', 'Fraction']]).transform('median')
        dataset['FractionMedian'] = dataset['RunMedian'].groupby(dataset['Fraction']).transform('median')
        dataset[NORM_INTENSITY] = dataset[INTENSITY] - dataset['RunMedian'] + dataset['FractionMedian']
        return dataset

    elif scaling_method == 'qnorm':
        # pivot to have one col per sample
        normalize_df = pd.pivot_table(dataset, index=[CONDITION, PROTEIN_NAME, PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, FRACTION, RUN, BIOREPLICATE, RT, SEARCH_ENGINE],
                                      columns=class_field, values=field, aggfunc={field: np.mean})
        if indisk:
            hdf_name = secrets.token_hex(nbytes=16) + ".parquet"
            normalize_df.to_parquet(hdf_name)
            df = pd.read_parquet(hdf_name)
            normalize_df = qnorm.quantile_normalize(df, ncpus=4)
        else:
            normalize_df = qnorm.quantile_normalize(normalize_df, axis=1)
        normalize_df = normalize_df.reset_index()
        normalize_df = normalize_df.melt(id_vars=[CONDITION, PROTEIN_NAME, PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, FRACTION, RUN, BIOREPLICATE, RT, SEARCH_ENGINE])
        normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
        print(dataset.head())
        return normalize_df

    #
    # # normalize the intensities across all samples
    # if class_field == "all":
    #     normalize_df = dataset[[field]]
    #     dataset[NORM_INTENSITY] = scaler.fit_transform(normalize_df)
    # else:
    #     # normalize taking into account samples
    #     normalize_df = dataset[[PEPTIDE_SEQUENCE, CONDITION, field, class_field]]
    #     # group peptide + charge into a single peptide intensity using the mean.
    #     normalize_df = pd.pivot_table(normalize_df, values=field, index=[PEPTIDE_SEQUENCE, CONDITION],
    #                                   columns=class_field,
    #                                   aggfunc={field: np.mean})
    #
    #     print(normalize_df.head())
    #
    #     # Remove all peptides in less than 30% of the samples.
    #     if remove_peptides:
    #         normalize_df = remove_missing_values(normalize_df=normalize_df, ratio = 0.3)
    #
    #     # Imputation of the values using KNNImputer
    #     # (https://scikit-learn.org/0.16/modules/generated/sklearn.preprocessing.Imputer.html)
    #     if imputation:
    #         imputer = SimpleImputer()
    #         if imputation_method != "simple":
    #             imputer = KNNImputer(n_neighbors=2, weights="uniform")
    #         normalized_matrix = imputer.fit_transform(normalize_df)
    #         if scaling_method != "quantile":
    #             normalized_matrix = scaler.fit_transform(normalized_matrix)
    #         else:
    #             normalized_matrix = qnorm.quantile_normalize(normalized_matrix)
    #     else:
    #         normalized_matrix = scaler.fit_transform(normalize_df)
    #
    #     normalize_df[:] = normalized_matrix
    #     normalize_df = normalize_df.reset_index()
    #     normalize_df = normalize_df.melt(id_vars=[PEPTIDE_SEQUENCE, CONDITION])
    #     normalize_df.rename(columns={'value': NORM_INTENSITY}, inplace=True)
    #     dataset = pd.merge(dataset, normalize_df, how='left', on=[PEPTIDE_SEQUENCE, CONDITION, class_field])
    #     dataset.pop(NORM_INTENSITY + "_x")
    #     dataset.rename(
    #         columns={NORM_INTENSITY + "_y": NORM_INTENSITY}, inplace=True)

    return dataset

def get_peptidoform_normalize_intensities(dataset: DataFrame, higher_intensity: bool = True) -> DataFrame:
    """
    Select the best peptidoform for the same sample and the same replicates. A peptidoform is the combination of
    a (PeptideSequence + Modifications) + Charge state.
    :param dataset: dataset including all properties
    :param higher_intensity: select based on normalize intensity, if false based on best scored peptide
    :return:
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    if higher_intensity:
        dataset = dataset.loc[dataset.groupby([PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, BIOREPLICATE])[NORM_INTENSITY].idxmax()].reset_index(drop=True)
    else:
        dataset = dataset.loc[dataset.groupby([PEPTIDE_SEQUENCE, PEPTIDE_CHARGE, SAMPLE_ID, BIOREPLICATE])[
            SEARCH_ENGINE].idxmax()].reset_index(drop=True)
    return dataset

def sum_peptidoform_intensities(dataset: DataFrame) -> DataFrame:
    """
    Sum the peptidoform intensities for all peptidofrom across replicates of the same sample.

    :param dataset: Dataframe to be analyzed
    :return: dataframe with the intensities
    """
    dataset = dataset[dataset[NORM_INTENSITY].notna()]
    dataset = dataset.groupby([CONDITION, PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID, BIOREPLICATE])[NORM_INTENSITY].sum()
    dataset = dataset.reset_index()
    return dataset

def average_peptide_intensities(dataset: DataFrame) -> DataFrame:
    """
    Median the intensities of all the peptidoforms for an specific peptide sample combination.
    :param dataset: Dataframe containing all the peptidoforms
    :return: New dataframe
    """
    dataset = dataset.groupby([CONDITION, PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID])[NORM_INTENSITY].median()
    dataset = dataset.reset_index()
    return dataset

def compute_correlations_cross_samples(dataset: DataFrame) -> None:
    normalize_df = pd.pivot_table(dataset,
                                  index=[CONDITION, PROTEIN_NAME, PEPTIDE_CANONICAL],
                                  columns=SAMPLE_ID, values=NORM_INTENSITY)
    print(normalize_df[list(normalize_df.columns)].corr())

def filter_by_peptide_length(row):
    peptide = row[PEPTIDE_SEQUENCE]
    canonical = get_canonical_peptide(peptide)
    if len(canonical) > 7:
        return True
    else:
        return False

@click.command()
@click.option("--peptides", help="Peptides files from the peptide file generation tool")
@click.option("--contaminants", help="Contaminants and high abundant proteins to be removed")
@click.option("--routliers", help="Remove outliers from the peptide table", is_flag=True)
@click.option("--output", help="Peptide intensity normalized")
@click.option('--nmethod', help="Normalization method used to normalize intensities for all samples (options: quantile, robusts, standard)", default="quantile")
@click.option("--compress", help="Read the input peptides file in compress gzip file", is_flag= True)
@click.option("--log2", help="Transform to log2 the peptide intensity values before normalization", is_flag=True)
@click.option("--violin", help="Use violin plot instead of boxplot for distribution representations", is_flag=True)
@click.option("--show", help="Show the plots for each steps", is_flag=True)
@click.option("--quality", help="Apply quality control to peptide information", is_flag=True)
@click.option("--indisk", help="Perform the quantile normalization in disk, not in memory", is_flag = True)
@click.option("--verbose",
              help="Print addition information about the distributions of the intensities, number of peptides remove after normalization, etc.",
              is_flag=True)
def peptide_normalization(peptides: str, contaminants: str, routliers: bool, output: str, nmethod: str, compress: bool, log2: bool,
                          violin: bool, show: bool, quality: bool, indisk: bool, verbose: bool) -> None:

    if peptides is None or output is None:
        print_help_msg(peptide_normalization)
        exit(1)

    output_file_prefix = os.path.splitext(output)[0]

    pd.set_option('display.max_columns', None)
    # TODO infer from filename
    print("Loading data..")
    compression_method = 'gzip' if compress else None
    if compress:
        dataset_df = pd.read_csv(peptides, sep="\t", compression=compression_method)
    else:
        dataset_df = pd.read_csv(peptides, sep="\t")
    print_dataset_size(dataset_df, "Number of peptides: ", verbose)

    print("Logarithmize if specified..")
    dataset_df[NORM_INTENSITY] = np.log2(dataset_df[INTENSITY]) if log2 else dataset_df[INTENSITY]

    # Print the distribution of the original peptide intensities from quantms analysis
    if verbose:
        plot_distributions(dataset_df, INTENSITY, SAMPLE_ID, log2=not log2, file_name=output_file_prefix+"-{}-{}.pdf".format("distribution","non-normalize"), show = show)
        plot_box_plot(dataset_df, INTENSITY, SAMPLE_ID, log2=not log2,
                      title="Original peptide intensity distribution (no normalization)", violin=violin, file_name=output_file_prefix+"-{}-{}.pdf".format("boxplot","non-normalize"), show=show)

    # Remove high abundant and contaminants proteins and the outliers
    if contaminants is not None:
        print("Remove contaminants...")
        dataset_df = remove_contaminants_decoys(dataset_df, contaminants=contaminants)
    print_dataset_size(dataset_df, "Peptides after contaminants removal: ", verbose)

    if quality:
        print("Removing small peptides...")
        mask = dataset_df.apply(filter_by_peptide_length, axis=1)
        dataset_df = dataset_df[mask]
        print(dataset_df.head())
    print_dataset_size(dataset_df, "Peptides after small peptides (< 7AA) removal: ", verbose)

    if verbose:
        plot_distributions(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=not log2, file_name=output_file_prefix+"-{}-{}.pdf".format("distribution","contaminant-remove"), show = show)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=not log2,
                      title="Peptide intensity distribution after contaminants removal", violin=violin, file_name=output_file_prefix+"-{}-{}.pdf".format("boxplot","contaminant-remove"), show = show)

    print("Normalize intensities.. ")
    dataset_df = intensity_normalization(dataset_df, field=NORM_INTENSITY, class_field=SAMPLE_ID, scaling_method=nmethod, indisk=indisk)

    if verbose:
        log_after_norm = nmethod == "msstats" or nmethod == "qnorm" or ((nmethod == "quantile" or nmethod == "robust") and not log2)
        plot_distributions(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm, file_name=output_file_prefix+"-{}-{}.pdf".format("distribution","normalized"), show = show)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm,
                      title="Peptide intensity distribution after imputation, normalization method: " + nmethod, violin=True, file_name=output_file_prefix+"-{}-{}.pdf".format("boxplot","normalized"), show = show)

    print("Select the best peptidoform across fractions...")
    print("Number of peptides before peptidofrom selection: " + str(len(dataset_df.index)))
    dataset_df = get_peptidoform_normalize_intensities(dataset_df)
    print("Number of peptides after peptidofrom selection: " + str(len(dataset_df.index)))

    # Add the peptide sequence canonical without the modifications
    print("Add Canonical peptides to the dataframe...")
    dataset_df[PEPTIDE_CANONICAL] = dataset_df[PEPTIDE_SEQUENCE].apply(lambda x: get_canonical_peptide(x))

    print("Sum all peptidoforms per Sample...")
    print("Number of peptides before sum selection: " + str(len(dataset_df.index)))
    dataset_df = sum_peptidoform_intensities(dataset_df)
    print("Number of peptides after sum: " + str(len(dataset_df.index)))

    print("Average all peptidoforms per Peptide/Sample...")
    print("Number of peptides before average: " + str(len(dataset_df.index)))
    dataset_df = average_peptide_intensities(dataset_df)
    print("Number of peptides after average: " + str(len(dataset_df.index)))

    if verbose:
        log_after_norm = nmethod == "msstats" or nmethod == "qnorm" or ((nmethod == "quantile" or nmethod == "robust") and not log2)
        plot_distributions(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm, file_name=output_file_prefix+"-{}-{}.pdf".format("distribution","final"), show = show)
        plot_box_plot(dataset_df, NORM_INTENSITY, SAMPLE_ID, log2=log_after_norm,
                      title="Peptide intensity distribution after imputation, normalization method: " + nmethod, violin=True, file_name=output_file_prefix+"-{}-{}.pdf".format("boxplot","final"), show = show)

    compute_correlations_cross_samples(dataset_df)

    dataset_df.to_csv(output, index=False, sep='\t')

if __name__ == '__main__':
    peptide_normalization()
