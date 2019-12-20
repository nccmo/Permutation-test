#!/usr/bin/env python2

import gzip
import sys
import os
import subprocess
import logging
import argparse
import re
import random
import array
import csv
import itertools

import pysam
import pandas as pd
import numpy as np

import permutation_test
import filtering
import utils
from gene_specific_analysis import *

logger = logging.getLogger(__name__)  # module logger

#trinucleotide context
trinucs = list(set(
        [t1+t2+t3
         for t1 in 'ACTG'
         for t2 in 'CT'
         for t3 in 'ACTG']
    ))
strand_trinucs = list(set(
        [t0+t1+t2+t3
         for t0 in '+-'
         for t1 in 'ACTG'
         for t2 in 'CT'
         for t3 in 'ACTG']
    ))

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] # , 'Y']

codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*',
               'Splice_Site': 'Splice_Site'}

re_refID =re.compile(r'\((N[MR]_\d+)\)')

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

# count silent mutations of the dataframe

def count_silent_mutations(opts):

    ############################
    # read in necessary files
    ############################

    # read in position to trinucleotide file
    logger.info('reading pos_to_class_dataframe...')
    pos_to_nuc_df = pd.read_table('db/merged_pos_to_context_class.txt', sep='\t', names=('pos', 'strand', 'trinucleotide', 'expression', 'replication_time', 'exp_group', 'rep_group', 'class'))
    pos_to_nuc_df['strand_trinucs'] = pos_to_nuc_df['strand'].str.cat(pos_to_nuc_df['trinucleotide'], sep='')

    del pos_to_nuc_df['strand']
    del pos_to_nuc_df['trinucleotide']
    del pos_to_nuc_df['expression']
    del pos_to_nuc_df['replication_time']
    del pos_to_nuc_df['exp_group']
    del pos_to_nuc_df['rep_group']
    
    for i,k in itertools.product(range(1,16), range(1,16)):
        tmp = pos_to_nuc_df['class'] == 'exp_class' + str(i) + '_rep_class' + str(k)
        exec('exp_class' + str(i) + '_rep_class' + str(k) + '_num = tmp')
        exec('tmp2 = str(exp_class' + str(i) + '_rep_class' + str(k) + '_num.sum())')
        print 'exp_class' + str(i) + '_rep_class' + str(k) + '\t' + tmp2
    
    #############################
    # modify mutation dataframe
    #############################

    # get mutation df
    mut_df = pd.read_csv(opts['maf_file'], sep='\t')
    orig_num_mut = len(mut_df)

    # rename columns to fit my internal column names
    rename_dict = {
        'Hugo_Symbol': 'Gene',
        'Tumor_Sample_Barcode': 'Tumor_Sample',
        'Tumor_Seq_Allele2' : 'Tumor_Allele',
        'Tumor_Seq_Allele' : 'Tumor_Allele',
    }
    mut_df.rename(columns=rename_dict, inplace=True)

    # drop rows with missing info
    na_cols = ['Gene', 'Tumor_Allele', 'Start_Position', 'Chromosome']
    mut_df = mut_df.dropna(subset=na_cols)
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Dropped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    #############################
    #. SNV dataframe
    #############################

    # select valid silent single nucleotide variants only and corrects for 1-based coordinates!! (important)
    silent_snv_df = filtering.silent_mutation_df(mut_df, opts['unique'])

    # get chromosome-position
    silent_snv_df['Chrom_Pos'] = silent_snv_df['Chromosome'] + ':' + silent_snv_df['Start_Position'].astype(str)

    # remove SNVs of non-coding regions
    orig_len = len(silent_snv_df) 
    removed_silent_snv_df = silent_snv_df[~silent_snv_df['Chrom_Pos'].isin(pos_to_nuc_df['pos'])]
    removed_silent_snv_df.to_csv("db/removed_silent_snv_df.txt", sep = '\t', header =False, index=False)
    silent_snv_df = silent_snv_df[silent_snv_df['Chrom_Pos'].isin(pos_to_nuc_df['pos'])]
    after_len = len(silent_snv_df)

    log_msg = ('Kept {num_left} coding silent mutations '
               'after dropping {num_dropped} non-coding mutations.'.format(num_left=after_len, num_dropped=orig_len-after_len))
    logger.info(log_msg)

    silent_snv_df_class = pd.merge(silent_snv_df, pos_to_nuc_df, left_on='Chrom_Pos', right_on='pos')

    for i,k in itertools.product(range(1,16), range(1,16)):
        tmp = silent_snv_df_class['class'] == 'exp_class' + str(i) + '_rep_class' + str(k)
        exec('exp_class' + str(i) + '_rep_class' + str(k) + ' = tmp')
        exec('tmp2 = str(exp_class' + str(i) + '_rep_class' + str(k) + '.sum())')
        print 'exp_class' + str(i) + '_rep_class' + str(k) + '\t' + tmp2

    print 'total' + '\t' + str(len(silent_snv_df_class))

    class_to_freq_dict = {}

    for i,k in itertools.product(range(1,16), range(1,16)):
        exec('tmp = exp_class' + str(i) + '_rep_class' + str(k) + '_num.sum() == 0')
        if tmp == True:
            print 'exp_class' + str(i) + '_rep_class' + str(k) + '\t' + str(0)
            class_to_freq_dict['exp_class' + str(i) + '_rep_class' + str(k)] = 0
        else:
            exec('tmp = float(1000.0 * exp_class' + str(i) + '_rep_class' + str(k) + '.sum() / exp_class' + str(i) + '_rep_class' + str(k) + '_num.sum())')
            print 'exp_class' + str(i) + '_rep_class' + str(k) + '\t' + str(tmp)
            class_to_freq_dict['exp_class' + str(i) + '_rep_class' + str(k)] = tmp

    pos_to_nuc_df['class_freq'] = pos_to_nuc_df['class'].apply(class_to_freq_dict.get)
    del pos_to_nuc_df['class']

    pos_to_nuc_df.to_csv("db/merged_pos_to_context_class_final.txt", sep = '\t', header =False, index=False)

    for k in strand_trinucs:
        tmp_df = pos_to_nuc_df[pos_to_nuc_df['strand_trinucs'] == k]
        tmp_df.to_csv("db/" + k + "_data.txt", sep = '\t', header =False, index=False)


def parse_arguments():
    # make a parser
    info = 'Count silent mutations'
    parent_parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parent_parser.add_argument('maf_file', default = None, type = str,
                               help = "the path to the maf file")

    help_str = ('Only keep unique mutations for each tumor sample. '
                'Mutations reported from heterogeneous sources may contain'
                ' duplicates, e.g. a tumor sample was sequenced twice.')

    parent_parser.add_argument('--unique', action='store_true', default=True,
                               help=help_str)

    args = parent_parser.parse_args()
    opts = vars(args)
    return opts


def cli_main():
    log_file = 'stdout'  # auto-name the log file
    log_level = ''
    utils.start_logging(log_file=log_file,
                        log_level='INFO',
                        verbose=False)  # start logging
    opts = parse_arguments()
    count_silent_mutations(opts)


if __name__ == "__main__":
    cli_main()