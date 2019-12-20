#!/usr/bin/env python2

import logging
import utils

logger = logging.getLogger(__name__)  # module logger

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

def snv_mutation_df(mutation_df, only_unique=True):

    # only keep allowed mutation types
    orig_len = len(mutation_df)  # number of mutations before filtering

    # select only SNVs
    flag_valid_nuc = (mutation_df['Reference_Allele'].apply(utils.check_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(utils.check_valid_nuc))
    mutation_df = mutation_df[flag_valid_nuc] 
    mutation_df = mutation_df[mutation_df['Tumor_Allele'].apply(lambda x: len(x)==1)]
    mutation_df = mutation_df[mutation_df['Reference_Allele'].apply(lambda x: len(x)==1)]
    valid_len = len(mutation_df)

    log_message = ('Dropped {num_dropped} mutations after only keeping SNVs'.format(num_dropped=orig_len-valid_len))
    logger.info(log_message)

    if only_unique:
        dup_cols = ['Tumor_Sample', 'Chromosome', 'Start_Position',
                    'End_Position', 'Reference_Allele', 'Tumor_Allele']
        mutation_df = mutation_df.drop_duplicates(subset=dup_cols)

        # log results of de-duplication
        dedup_len = len(mutation_df)
        log_message = ('Dropped {num_dropped} mutations when removing '
                       'duplicates'.format(num_dropped=valid_len-dedup_len))
        logger.info(log_message)

    # correct for 1-based coordinates
    mutation_df['Start_Position'] = mutation_df['Start_Position'].astype(int) - 1
    return mutation_df



def other_mutation_df(mutation_df, only_unique=True):

    # only keep allowed mutation types
    orig_len = len(mutation_df)  # number of mutations before filtering

    # check if mutations are valid insertion
    flag_valid_nuc = (mutation_df['Reference_Allele'].apply(utils.check_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(utils.check_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(lambda x: len(x)==1) & \
                      mutation_df['Reference_Allele'].apply(lambda x: len(x)==1))
    #non-snv
    mutation_df = mutation_df[~flag_valid_nuc]
    valid_len = len(mutation_df)

    log_message = ('Dropped {num_dropped} SNV mutations'.format(num_dropped=orig_len-valid_len))
    logger.info(log_message)

    if only_unique:
        dup_cols = ['Tumor_Sample', 'Chromosome', 'Start_Position',
                    'End_Position', 'Reference_Allele', 'Tumor_Allele']
        mutation_df = mutation_df.drop_duplicates(subset=dup_cols)

        # log results of de-duplication
        dedup_len = len(mutation_df)
        log_message = ('Dropped {num_dropped} mutations when removing '
                   'duplicates'.format(num_dropped=valid_len-dedup_len))
        logger.info(log_message)

    # correct for 1-based coordinates
    mutation_df['Start_Position'] = mutation_df['Start_Position'].astype(int) - 1
    return mutation_df



def silent_mutation_df(mutation_df, only_unique=True):

    mutation_df = mutation_df[mutation_df.Variant_Classification.isin(['Silent'])]  # only keep silent SNV

    # log the number of dropped mutations
    log_message = ('{num_dropped} silent mutations were detected in the dataset.'.format(num_dropped=len(mutation_df)))
    logger.info(log_message)

    # select only SNVs
    orig_len = len(mutation_df)
    flag_valid_nuc = (mutation_df['Reference_Allele'].apply(utils.check_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(utils.check_valid_nuc))
    mutation_df = mutation_df[flag_valid_nuc]
    mutation_df = mutation_df[mutation_df['Tumor_Allele'].apply(lambda x: len(x)==1)]
    mutation_df = mutation_df[mutation_df['Reference_Allele'].apply(lambda x: len(x)==1)]
    valid_len = len(mutation_df)

    log_message = ('Dropped {num_dropped} mutations after only keeping SNVs'.format(num_dropped=orig_len-valid_len))
    logger.info(log_message)

    if only_unique:
        dup_cols = ['Tumor_Sample', 'Chromosome', 'Start_Position',
                    'End_Position', 'Reference_Allele', 'Tumor_Allele']
        mutation_df = mutation_df.drop_duplicates(subset=dup_cols)

        # log results of de-duplication
        dedup_len = len(mutation_df)
        log_message = ('Dropped {num_dropped} mutations when removing '
                   'duplicates'.format(num_dropped=valid_len-dedup_len))
        logger.info(log_message)

    # correct for 1-based coordinates
    mutation_df['Start_Position'] = mutation_df['Start_Position'].astype(int) - 1
    return mutation_df
