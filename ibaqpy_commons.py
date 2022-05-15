from pandas import DataFrame

PROTEIN_NAME = 'ProteinName'
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
