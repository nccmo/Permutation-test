#! /usr/bin/env python2

import logging
import pandas as pd
import sys
logger = logging.getLogger(__name__)  # module logger

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

def gene_analysis(mut_df, interested_gene, pos_to_nuc_keys, maf_file):

    # select samples with non-synonymous mutation / INDELs in the coding region of the gene
    gene_mut_df = mut_df[mut_df['Gene'].isin([interested_gene])]
    gene_mut_df = gene_mut_df[~gene_mut_df['Variant_Classification'].isin(['Silent'])]
    
    gene_mut_df['Start_Position'] = gene_mut_df['Start_Position'].astype(int) - 1
    gene_mut_df['End_Position'] = gene_mut_df['End_Position'].astype(int) - 1

    gene_mut_df['Chrom_Pos'] = gene_mut_df['Chromosome'] + ':' + gene_mut_df['Start_Position'].astype(str)
    gene_mut_df['Chrom_Pos2'] = gene_mut_df['Chromosome'] + ':' + gene_mut_df['End_Position'].astype(str)

    gene_mut_df = gene_mut_df[(gene_mut_df['Chrom_Pos'].isin(pos_to_nuc_keys)) | (gene_mut_df['Chrom_Pos2'].isin(pos_to_nuc_keys))]

    gene_mut_samples = gene_mut_df['Tumor_Sample'].unique()

    log_msg = ('{sample_number} samples had non-synonymous mutations and/or INDELs in '
               '{target_gene} (dataset: {maf_file})'.format(sample_number=len(gene_mut_samples), target_gene=interested_gene, maf_file=maf_file))
    logger.info(log_msg)

    gene_multiple_mut_samples = gene_mut_df['Tumor_Sample'].value_counts()[gene_mut_df['Tumor_Sample'].value_counts()>=2]
    log_msg = ('{sample_number} samples had multiple mutations and/or INDELs in '
               '{target_gene} (dataset: {maf_file})'.format(sample_number=len(gene_multiple_mut_samples), target_gene=interested_gene, maf_file=maf_file))
    logger.info(log_msg)

    # mutation df of the gene_mut_samples
    analysis_df = mut_df[mut_df['Tumor_Sample'].isin(gene_mut_samples)]
    
    # remove interested_gene_mutations
    # inluding silents, UTRs.. etc
    analysis_df = analysis_df[~analysis_df['Gene'].isin([interested_gene])]

    return(analysis_df)