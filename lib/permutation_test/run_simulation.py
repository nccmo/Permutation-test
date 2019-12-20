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

import pysam
import pandas as pd
import numpy as np

import permutation_test
import filtering
import utils
from gene_specific_analysis import *
from collections import OrderedDict

logger = logging.getLogger(__name__)  # module logger

strand_trinucs = list(set(
        [t0+t1+t2+t3
         for t0 in '+-'
         for t1 in 'ACTG'
         for t2 in 'CT'
         for t3 in 'ACTG']
    ))

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X' ]# , 'Y']

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

base_pairing = {'A': 'T',
                'T': 'A',
                'a': 't',
                't': 'a',
                'C': 'G',
                'G': 'C',
                'c': 'g',
                'g': 'c',
                '-': '-',  # some people denote indels with '-'
                'n': 'n',
                'N': 'N'}

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

def weighted_choice(values, cdf, size=1):
    uniform_samples = np.random.sample(size)
    idx = cdf.searchsorted(uniform_samples, side='right')
    sample = values[idx]
    return sample[0]


def SNV_main(args, mut_df=None, frameshift_df=None):

    opts = vars(args)

    ############################
    # read in necessary files
    ############################

    # read in position to trinucleotide file
    logger.info('reading pos_to_nuc_dictionary...')
    # read in data frame
    pos_to_nuc_df = pd.read_table('db/merged_pos_to_context_class_final.txt', sep='\t', names=('pos', 'trinucleotide', 'coefficient'))

    if len(pos_to_nuc_df[pos_to_nuc_df['trinucleotide'].astype(str).str.len() != 4]) != 0:
        logger.info('something is wrong with reading pos_to_nuc_dictionary...')
        sys.exit()

    # make dictionary
    pos_to_nuc = pos_to_nuc_df.set_index('pos')['trinucleotide'].to_dict()
    pos_to_nuc_keys = pos_to_nuc.keys()

    # read in trinucleotide to position file
    logger.info('reading nuc_to_pos_dictionary...')
    nuc_to_pos_dict = {}
    nuc_to_cumsum_dict = {}
    for k in strand_trinucs:
        df = pd.read_table('db/' + k + '_data.txt', sep='\t', names=('pos', 'trinucleotide', 'coefficient'))
        nuc_to_pos_dict[k] = df['pos'].values

        p = df['coefficient'].values
        cdf = np.cumsum(p)
        cdf /= cdf[-1]
        nuc_to_cumsum_dict[k] = cdf

    # read in cancer gene file
    logger.info('reading cancer_gene_dictionary...')
    pos_to_codon_dict = {}
    pos_to_gene_dict = {}
    tmp_Chr_pos = ''
    with open('db/cancergene_pos_list_final.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split(' ')
            if len(F) < 4:
                continue
            if opts['gene'] and F[0] != opts['gene']:
                continue
            if len(F[3]) != 3 | len(F[3]) != 11:
                logger.info("codon frame error...: {0}, {1}, {2}".format(F[0], F[1], F[3]))
                sys.exit()

            pos_to_gene_dict.setdefault(F[1], []).append(F[0])
            pos_to_codon_dict.setdefault(F[1], []).append(';'.join([F[2], F[3], F[4]]))

    hin.close()
    pos_to_gene_dict_keys = pos_to_gene_dict.keys()


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
    na_cols = ['Gene', 'Reference_Allele', 'Tumor_Allele', 'Start_Position', 'Chromosome']
    mut_df = mut_df.dropna(subset=na_cols)
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    if opts['gene']:
        mut_df = gene_analysis(mut_df, opts['gene'], pos_to_nuc_keys, opts['maf_file'])


    #############################
    #. SNV dataframe
    #############################

    # select valid single nucleotide variants only and corrects for 1-based coordinates!! (important)
    snv_df = filtering.snv_mutation_df(mut_df, opts['unique'])
    # get chromosome-position
    snv_df['Chrom_Pos'] = snv_df['Chromosome'] + ':' + snv_df['Start_Position'].astype(str)

    # remove SNVs of non-coding regions
    orig_len = len(snv_df) 
    snv_df = snv_df[snv_df['Chrom_Pos'].isin(pos_to_nuc_keys)]
    after_len = len(snv_df)

    log_msg = ('Dropped {num_dropped} non-coding SNV mutations.'.format(num_dropped=orig_len-after_len))
    logger.info(log_msg)


    #############################
    #. SNV check
    ############################

    # get trincleotide context
    snv_df['trinucleotide'] =snv_df['Chrom_Pos'].apply(lambda x:pos_to_nuc[x])
    snv_df['trinucleotide'] = snv_df['trinucleotide'].astype('category')
    snv_df['Chrom_Pos'] = snv_df['Chrom_Pos'].astype('category')

    # check if the mutation is in the gene_list
    snv_df['gene'] = snv_df['Chrom_Pos'].map(pos_to_gene_dict)

    tmp_snv_df = snv_df.dropna(subset=['gene'])
    outcome = []
    chr_pos = tmp_snv_df.Chrom_Pos.values
    t_allele = tmp_snv_df.Tumor_Allele.values
    n_allele = tmp_snv_df.Reference_Allele.values

    # check if the mutation is synonymous / non-synonymous / splice site
    for idx in range(tmp_snv_df.shape[0]):
        tmp_outcome = []
        # there are genes with different reading frames
        for item in pos_to_codon_dict[chr_pos[idx]]:
            
            pos_in_codon = item.split(';')[0]
            codon_seq = item.split(';')[1]
            strand = item.split(';')[2]

            if pos_in_codon == 'splice_site':
                tmp_outcome.append('splice_site')
                continue
            
            # check if base change causes amino acid change
            codon_seq_list = list(codon_seq)

            if codon_seq_list[int(pos_in_codon)] == n_allele[idx]:
                codon_seq_list[int(pos_in_codon)] = t_allele[idx]
            elif codon_seq_list[int(pos_in_codon)] == utils.reverse_complement(n_allele[idx]):
                codon_seq_list[int(pos_in_codon)] = utils.reverse_complement(t_allele[idx])
            else:
                print "error: " + chr_pos[idx] + pos_to_codon_dict[chr_pos[idx]]

            new_codon_seq = ''.join(codon_seq_list)

            if codon_table[codon_seq] == codon_table[new_codon_seq]:
                tmp_outcome.append('synonymous')
            else:
                tmp_outcome.append('non-synonymous')

        outcome.append(':'.join(tmp_outcome))

    tmp_snv_df['outcome'] = outcome
    tmp_snv_df['original'] = 'original'
    tmp_snv_df.to_csv(opts['output_prefix'] + '.final_snv_result.csv', columns=['Tumor_Sample', 'original', 'Gene', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele', 'Chrom:Pos', 'gene', 'outcome'], index=False)


    #############################
    #. SNV simulation
    ############################

    max_num_sim = opts['simulation_number'] # number of simulations

    for num_sim in range(max_num_sim):

        log_msg = ('Performing simulation {num_simulation}...'.format(num_simulation=num_sim+1))
        logger.info(log_msg)

        # randomization
        trinuc = snv_df['trinucleotide'].values
        new_pos_list = []
        for idx in range(snv_df.shape[0]):
            wr = weighted_choice(nuc_to_pos_dict[trinuc[idx]], nuc_to_cumsum_dict[trinuc[idx]])
            new_pos_list.append(wr)
            
        snv_df['New_chr_pos'] = new_pos_list

        # check if new chr_pos is in gene_list
        snv_df['New_gene'] = snv_df['New_chr_pos'].map(pos_to_gene_dict)

        tmp_snv_df = snv_df.dropna(subset=['New_gene'])
        #print tmp_snv_df

        outcome = []
        chr_pos = tmp_snv_df.New_chr_pos.values
        t_allele = tmp_snv_df.Tumor_Allele.values
        n_allele = tmp_snv_df.Reference_Allele.values

        for idx in range(tmp_snv_df.shape[0]):
            tmp_outcome = []
            for item in pos_to_codon_dict[chr_pos[idx]]:
                
                pos_in_codon = item.split(';')[0]
                codon_seq = item.split(';')[1]
                strand = item.split(';')[2]

                if pos_in_codon == 'splice_site':
                    tmp_outcome.append('splice_site')
                    continue

                codon_seq_list = list(codon_seq)

                if codon_seq_list[int(pos_in_codon)] == n_allele[idx]:
                    codon_seq_list[int(pos_in_codon)] = t_allele[idx]
                elif codon_seq_list[int(pos_in_codon)] == utils.reverse_complement(n_allele[idx]):
                    codon_seq_list[int(pos_in_codon)] = utils.reverse_complement(t_allele[idx])
                else:
                    print "error: " + chr_pos[idx] + pos_to_codon_dict[chr_pos[idx]]

                new_codon_seq = ''.join(codon_seq_list)

                if codon_table[codon_seq] == codon_table[new_codon_seq]:
                    tmp_outcome.append('synonymous')
                else:
                    tmp_outcome.append('non-synonymous')

            outcome.append(':'.join(tmp_outcome))
                
        tmp_snv_df['New_outcome'] = outcome
        tmp_snv_df['simulation_num'] = 'simulation' + str(int(num_sim)+1)
        tmp_snv_df.to_csv(opts['output_prefix'] + '.final_snv_result.csv', columns=['Tumor_Sample', 'simulation_num', 'Gene', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele', 'New_chr_pos', 'New_gene', 'New_outcome'], mode='a', header=False, index=False)

    log_msg = ('Successfully finished. gene:{gene}, maf:{maf}'.format(gene=opts['gene'], maf=opts['maf_file']))
    logger.info(log_msg)


def non_SNV_main(args):
    
    opts = vars(args)

    ############################
    # read in necessary files
    ############################

    # read in position to trinucleotide file
    logger.info('reading pos_to_nuc_dictionary...')
    # read in data frame
    pos_to_nuc_df = pd.read_table('db/merged_pos_to_context_class_final.txt', sep='\t', names=('pos', 'trinucleotide', 'coefficient'))
    # make dictionary
    pos_to_nuc = pos_to_nuc_df.set_index('pos')['trinucleotide'].to_dict()
    pos_to_nuc_keys = pos_to_nuc.keys()

    # read in cancer gene file
    logger.info('reading cancer_gene_dictionary...')
    pos_to_codon_dict = {}
    pos_to_gene_dict = {}
    tmp_Chr_pos = ''
    with open('db/cancergene_pos_list_final.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split(' ')
            if len(F) < 4:
                continue
            if opts['gene'] and F[0] != opts['gene']:
                continue

            pos_to_gene_dict.setdefault(F[1], []).append(F[0])
            pos_to_codon_dict.setdefault(F[1], []).append(';'.join([F[2], F[3], F[4]]))
    hin.close()
    pos_to_gene_dict_keys = pos_to_gene_dict.keys()


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
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    if opts['gene']:
        mut_df = gene_analysis(mut_df, opts['gene'], pos_to_nuc_keys, opts['maf_file'])


    #############################
    #. others dataframe
    #############################

    # others SNV
    others_df = filtering.other_mutation_df(mut_df, opts['unique'])
    # get chromosome-position
    others_df['Chrom_Pos'] = others_df['Chromosome'] + ':' + others_df['Start_Position'].astype(str)

    # remove non-SNVs of non-exonic regions
    orig_len = len(others_df) 
    others_df = others_df[others_df['Chrom_Pos'].isin(pos_to_nuc_keys)]
    after_len = len(others_df)

    log_msg = ('Dropped {num_dropped} non-coding non-SNV mutations.'.format(num_dropped=orig_len-after_len))
    logger.info(log_msg)

    others_df['Chrom_Pos'] = others_df['Chrom_Pos'].astype('category')


    #############################
    #. others check
    ############################
    others_df['gene'] = others_df['Chrom_Pos'].map(pos_to_gene_dict)
    tmp_others_df = others_df.dropna(subset=['gene'])

    tmp_others_df['outcome'] = tmp_others_df['Chrom_Pos'].apply(lambda x:'INDELs' if x in pos_to_gene_dict_keys else 'NA')
    tmp_others_df['original'] = 'original'
    tmp_others_df.to_csv(opts['output_prefix'] + '.final_non_snv_result.csv', columns=['Tumor_Sample', 'original', 'Gene', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele', 'Chrom:Pos', 'gene', 'outcome'], index=False)


    #############################
    #. others simulation
    #############################

    max_num_sim = opts['simulation_number'] # number of simulations

    for num_sim in range(max_num_sim):

        log_msg = ('Performing {num_simulation} simulation ...'.format(num_simulation=num_sim+1))
        logger.info(log_msg)

        # randomization
        others_df['New_chr_pos'] = [random.choice(pos_to_nuc_keys) for _ in range(others_df.shape[0])]

        # check gene
        others_df['New_gene'] = others_df['New_chr_pos'].map(pos_to_gene_dict)
        tmp_others_df = others_df.dropna(subset=['New_gene'])

        # check outcome
        tmp_others_df['New_outcome'] = tmp_others_df['New_chr_pos'].apply(lambda x:'INDELs' if x in pos_to_gene_dict_keys else 'NA')
        tmp_others_df['simulation_num'] = 'simulation' + str(int(num_sim)+1)

        tmp_others_df.to_csv(opts['output_prefix'] + '.final_non_snv_result.csv', columns=['Tumor_Sample', 'simulation_num', 'Gene', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Allele', 'New_chr_pos', 'New_gene', 'New_outcome'], mode='a', header=False, index=False)

    #others_df.to_csv(opts['output_prefix'] + '.non_snv_df.csv')
    log_msg = ('Successfully finished. gene:{gene}, maf:{maf}'.format(gene=opts['gene'], maf=opts['maf_file']))
    logger.info(log_msg)



def parse_arguments():

    # make a parser
    info = 'script for permutation_test'
    parser = argparse.ArgumentParser(prog = "permutation_test", description=info, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    ####################

    # SNVs
    SNV = subparsers.add_parser("SNV", 
                                  help = "permutation test on SNVs",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # necessary input file
    SNV.add_argument("maf_file", default = None, type = str,
                        help = "the path to the maf file")

    SNV.add_argument("output_prefix", metavar = "output_prefix", default = None, type = str, 
                        help = "the prefix of the output")

    help_str = ('Only keep unique mutations for each tumor sample. '
                'Mutations reported from heterogeneous sources may contain'
                ' duplicates, e.g. a tumor sample was sequenced twice.')

    SNV.add_argument('--unique', action='store_true', default=True,
                     help=help_str)

    help_str = ('this command is for gene_specific analysis')

    SNV.add_argument('-g', '--gene', default=None, type = str,
                     help=help_str)

    help_str = ('Number of Monte Carlos simulations. default is 100.')

    SNV.add_argument('-n', '--simulation_number', default = 100, type = int,
                     help = help_str)

    SNV.add_argument('-f', '--fasta_file', default = "/home/w3varann/genomon_pipeline-2.4.0/database/GRCh37/GRCh37.fa", type = str,
                        help = "the path to the referece fasta file")

    SNV.set_defaults(func = SNV_main)

    ####################

    # non_SNVs
    non_SNV = subparsers.add_parser("non_SNV", 
                                  help = "permutation test on non_SNVs",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # necessary input file
    non_SNV.add_argument("maf_file", default = None, type = str,
                        help = "the path to the maf file")

    non_SNV.add_argument("output_prefix", metavar = "output_prefix", default = None, type = str, 
                        help = "the prefix of the output")

    help_str = ('Only keep unique mutations for each tumor sample. '
                'Mutations reported from heterogeneous sources may contain'
                ' duplicates, e.g. a tumor sample was sequenced twice.')

    non_SNV.add_argument('--unique', action='store_true', default=True,
                         help=help_str)

    help_str = ('this command is for gene_specific analysis')

    non_SNV.add_argument('-g', '--gene', default=None, type = str,
                        help=help_str)

    help_str = ('Number of Monte Carlos simulations. default is 100.')

    non_SNV.add_argument('-n', '--simulation_number', default = 100, type = int,
                         help = help_str)

    # logging arguments
    non_SNV.add_argument('-f', '--fasta_file', default = "/home/w3varann/genomon_pipeline-2.4.0/database/GRCh37/GRCh37.fa", type = str,
                        help = "the path to the referece fasta file")

    non_SNV.set_defaults(func = non_SNV_main)

    args = parser.parse_args()
    args.func(args)


def cli_main():
    log_file = 'stdout'  # auto-name the log file
    log_level = ''
    utils.start_logging(log_file=log_file,
                        log_level='INFO',
                        verbose=False)  # start logging
    parse_arguments()

if __name__ == "__main__":
    cli_main()
