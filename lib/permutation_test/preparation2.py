#!/usr/bin/env python2

import gzip
import re
import pysam
import sys
import pandas as pd
import subprocess
import logging
import utils

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] # , 'Y']
re_refID =re.compile(r'\((N[MR]_\d+)\)')
fa = pysam.Fastafile('/home/w3varann/genomon_pipeline-2.4.0/database/GRCh37/GRCh37.fa')
trinucs = list(set(
        [t1+t2+t3
         for t1 in 'ACTG'
         for t2 in 'CT'
         for t3 in 'ACTG']
    ))
strand_pairing = {'+': '-', '-': '+'}
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

logger = logging.getLogger(__name__)  # module logger

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

def gencode_codon_list(target_gene, output):

    gencode_df = pd.read_table('db/gencode_coding.modified.bed', names=('chr', 'start', 'end', 'ID', 'type', 'strand', 'gene', 'order', 'sum'))
    gencode_df = gencode_df[gencode_df['gene'] == target_gene]
    gencode_df = gencode_df[gencode_df['chr'] != 'chrY']
    ID_list = gencode_df['ID'].unique()

    # select only unique ID
    ID_list = list(set(ID_list))

    # get sequence for each refID
    hin = open("db/gencode_coding.modified.bed", 'r')
    gene_seq_dict = {}
    for line in hin:
        F = line.rstrip('\n').split('\t')
        chrom = F[0].replace('chr','')
        if chrom not in chroms:
            continue
        if F[4] != "coding":
            continue
        coding_start = int(F[1])
        coding_end   = int(F[2])
        strand = F[5]
        ID = F[3]
        gene = F[6]
        if gene != target_gene:
            continue
        for item in ID_list:
            if item not in gene_seq_dict: gene_seq_dict[item] = ''
            if ID == item:
                exon_seq = fa.fetch(reference=chrom, start=coding_start, end =coding_end).upper()
                gene_seq_dict[item] = gene_seq_dict[item] + exon_seq
    hin.close()

    hout = open(output, 'a')
    # get codon information from refgene
    hin = open("db/gencode_coding.modified.bed", 'r')
    exon_length_dict = {}
    for line in hin:
        F = line.rstrip('\n').split('\t')
        chrom = F[0].replace('chr','')
        if chrom not in chroms:
            continue
        if F[4] == "coding":
            coding_start = int(F[1])
            coding_end   = int(F[2])
            strand = F[5]
            ID = F[3]
            gene = F[6]
            if gene != target_gene:
                continue
            for item in ID_list:
                if item != ID: continue
                if item not in exon_length_dict: exon_length_dict[item] = 0
                if strand == '+':
                    for pos in range(coding_end - coding_start):
                        # relative pos in the gene
                        pos2 = pos + exon_length_dict[item]
                        codon_pos = pos2 // 3
                        codon_start = codon_pos * 3
                        pos_in_codon = pos2 % 3
                        print >> hout, gene, item, ':'.join([chrom, str(coding_start + pos)]), pos_in_codon, gene_seq_dict[item][codon_start:(codon_start+3)], strand
                    exon_length_dict[item] = exon_length_dict[item] + (coding_end - coding_start)

                if strand == '-':
                    for pos in range(coding_end - coding_start):
                        # relative pos in the gene
                        pos2 = pos + exon_length_dict[item]
                        codon_pos = pos2 // 3
                        codon_start = codon_pos * 3
                        pos_in_codon = pos2 % 3
                        print >> hout, gene, item, ':'.join([chrom, str(coding_start + pos)]), 2 - (pos_in_codon), utils.reverse_complement(gene_seq_dict[item][codon_start:(codon_start+3)]), strand
                    exon_length_dict[item] = exon_length_dict[item] + (coding_end - coding_start)

        elif F[4] == "intron":
            start = int(F[1])
            end   = int(F[2])
            strand = F[5]
            ID = F[3]
            gene = F[6]
            if gene != target_gene:
                continue
            for pos in (start, start + 1, end - 2, end - 1):
                print >> hout, gene, ID, ':'.join([chrom, str(pos)]), 'splice_site', 'splice_site', strand
    hin.close()
    hout.close()


####################################

def main():
    with open('input_db/all_cancergene_list.txt', 'r') as hin:
       for line in hin:
            F = line.rstrip('\n')
            if F.startswith('#'):
                continue
            logger.info('creating {F} gene infomation ...'.format(F=F))
            gencode_codon_list(F, 'db/cancergene_pos_list.txt')
    hin.close()

    df = pd.read_table("db/cancergene_pos_list.txt", sep='\s+', header=None)
    incomp_IDs = df[(df[4].astype(str).str.len() == 1) | (df[4].astype(str).str.len() == 2)][1].unique()
    df = df[~df[1].isin(incomp_IDs)]

    del df[1]
    df = df.drop_duplicates()
    df.to_csv('db/cancergene_pos_list_filtered.txt', sep = ' ', header =False, index=False)

    # sort
    hout = open('db/cancergene_pos_list_final.txt', 'w')
    s_ret = subprocess.call(["sort", "-k1", 'db/cancergene_pos_list_filtered.txt'], stdout = hout)

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting merged junction file"
        sys.exit(1)

    hin.close()
    hout.close()

if __name__ == "__main__":

    log_file = 'stdout'  # auto-name the log file
    log_level = ''
    utils.start_logging(log_file=log_file,
                        log_level='INFO',
                        verbose=False)  # start logging
    main()