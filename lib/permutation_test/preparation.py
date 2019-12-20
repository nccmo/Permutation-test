#!/usr/bin/env python2

import annot_utils
import pysam
import os
import gzip
import subprocess
import pandas as pd
import numpy as np
import sys
import re
import math
import logging
import utils

logger = logging.getLogger(__name__)  # module logger

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] #, 'Y']
trinucs = list(set(
        [t1+t2+t3
         for t1 in 'ACTG'
         for t2 in 'CT'
         for t3 in 'ACTG']
    ))
strand_pairing = {'+': '-', '-': '+'}

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

def db_preparation_main():

    #start logging
    log_file = 'stdout'  # auto-name the log file
    log_level = ''
    utils.start_logging(log_file=log_file,
                        log_level='INFO',
                        verbose=False)  # start logging

    # change directory to db/
    if not os.path.exists('db/'):
        os.mkdir('db/')
    os.chdir('db/')

    #####
    # make necessary bed files using annot_utils
    #####

    logger.info('making annot_utils files...')

    # make refseq gene coding bed file (refseq)
    annot_utils.coding.make_coding_info("gencode_coding.bed.gz", "gencode", "hg19", False, True)
    annot_utils.gene.make_gene_info("gencode_gene.bed.gz", "gencode", "hg19", False, True)

    logger.info('finished making annot_utils files...')

    # get trinucleotide context for each bed line
    # hack to index the FASTA file
    fa = pysam.Fastafile('/home/w3varann/genomon_pipeline-2.4.0/database/GRCh37/GRCh37.fa')

    # gencode dictionary
    
    hin = gzip.open('/home/ysaito/bin/annot_utils-0.2.1/annot_utils/data/hg19/wgEncodeGencodeCompV19.txt.gz', 'r')
    ID_dict = {}

    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene = str(F[12])
        genID = F[1]
        ID_dict[genID] = str(gene)
    hin.close()

    ######
    # Average Expression Data
    ######

    logger.info('reading in expression data...')
    
    expr_dict = {}

    with open("../input_db/CCLE_avg_exp.txt") as handle:
        # skip comment line(s)
        comment_line = handle.readline()
        header = handle.readline()
        for line in handle:
            F = line.rstrip('\n').split('\t')
            # consider genes with multiple transcripts
            if F[0] in expr_dict.keys():
                expr_dict[F[0]] = str(F[1]) + ":" + str(expr_dict[F[0]])
            else:
                expr_dict[F[0]] = F[1]
    handle.close()

    # get average expression for genes with multiple transcripts
    for k, v in expr_dict.items():
        if ":" in v:
            v_list = [float(n) for n in v.split(':')]
            expr_dict[k] = sum(v_list) / len(v_list)

    # get median expression of all genes
    all_median_expr = np.nanmedian([float(v) for v in expr_dict.values()])
    
    logger.info('finished reading in expression data')


    ######
    # Average Replication Time Data (for gene)
    ######

    logger.info('reading in replication_time data...')

    # skip comment line(s)
    with open("../input_db/replication_time.txt") as hIN:
        first_line = next(hIN)
        if first_line.startswith('#'):
            skip_rows = 1 if first_line.startswith('#') else 0
    hIN.close()

    # read in data frame
    rep_df = pd.read_table("../input_db/replication_time.txt", sep='\t', skiprows=skip_rows)
    # make dictionary
    rep_dict = rep_df.set_index('gene')['replication_time'].to_dict()


    ######
    # Average Replication Time Data (for position)
    ######

    # skip comment line(s)
    with open("../input_db/replication_time_position.txt") as hIN:
        first_line = next(hIN)
        if first_line.startswith('#'):
            skip_rows = 1 if first_line.startswith('#') else 0
    hIN.close()

    # read in data frame
    rep_pos_df = pd.read_table("../input_db/replication_time_position.txt", sep='\t', skiprows=skip_rows)

    rep_pos_df['chr_start'] = rep_pos_df['chr'].astype(str) + ':' + rep_pos_df['start'].astype(str)
    # make dictionary
    rep_pos_dict = rep_pos_df.set_index('chr_start')['replication_time'].to_dict()

    # get averaage replication time of all position
    all_median_rep_time = np.nanmedian([float(v) for v in rep_pos_dict.values()])

    logger.info('finished reading in replication_time data')


    ######
    # Fill in missing replication time data and expression data 
    ######

    logger.info('filling in missing replication time data and expression data...')

    with gzip.open('gencode_gene.bed.gz') as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            chrom = F[0].replace('chr','')
            if chrom not in chroms:
                continue
            gene_start = int(F[1])
            gene_end   = int(F[2])
            gene_mid_pos = (gene_start + gene_end) / 2
            gene = ID_dict[F[3]]

            if gene not in rep_dict.keys():
                gene_mid_pos_floor = int(math.floor(gene_mid_pos/100000)*100000 + 1)
                pos = chrom  + ':' + str(gene_mid_pos_floor)
                if pos in rep_pos_dict.keys():
                    rep_dict[gene] = rep_pos_dict[chrom + ':' + str(gene_mid_pos_floor)]
                else:
                    rep_dict[gene] = all_median_rep_time

            if gene not in expr_dict.keys():
                expr_dict[gene] = all_median_expr
    hIN.close()

    logger.info('finished filling in missing replication time data and expression data')


    ######
    # get trinucleotide context
    ######

    logger.info('getting trinucleotide context...')

    gencode_df = pd.read_table("gencode_coding.bed.gz", names=('chr', 'start', 'end', 'ID', 'type', 'strand'))
    # select only coding
    gencode_df = gencode_df[(gencode_df['type'] == 'coding') | (gencode_df['type'] == 'intron')]

    gencode_df['gene'] = gencode_df['ID'].map(ID_dict)

    # sort by ID
    gencode_df = gencode_df.sort_index(by = ['ID', 'chr', 'start', 'end'])
    #remove chrY
    gencode_df = gencode_df[gencode_df['chr'] != 'chrY']

    gencode_df.to_csv("gencode_coding.modified.bed", sep='\t', header=False, index=False)



    pos_to_context = open("merged_pos_to_context.txt", 'w')
    with open("gencode_coding.modified.bed", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chrom = F[0].replace('chr','')
            if chrom not in chroms:
                continue

            start = int(F[1])
            end   = int(F[2])
            length = end - start
            ID = F[3]
            gene = F[6]

            if F[4] == "coding":
                for i in range(length):
                    strand = F[5]
                    pos = start + i
                    try:
                        # 0-based
                        nucs = fa.fetch(reference=chrom,
                                        start=pos-1,
                                        end=pos+2).upper()
                        if nucs[1] == 'G':
                            strand = strand_pairing[strand]
                            nucs = utils.rev_seq(nucs)
                        if nucs[1] == 'A':
                            strand = strand_pairing[strand]
                            nucs = utils.rev_seq(nucs)

                    except Exception as inst:
                        logger.debug("{0}: {1}".format(type(inst), inst.args))
                    
                    if 'N' not in nucs:
                        print >> pos_to_context, '\t'.join([gene, chrom + ':' + str(pos), strand, nucs, str(expr_dict[gene]), str(rep_dict[gene])])

            if F[4] == "intron":
                for pos in (start, start + 1, end -2, end -1):
                    strand = F[5]
                    try:
                        # 0-based
                        nucs = fa.fetch(reference=chrom,
                                        start=pos-1,
                                        end=pos+2).upper()
                        if nucs[1] == 'G':
                            strand = strand_pairing[strand]
                            nucs = utils.rev_seq(nucs)
                        if nucs[1] == 'A':
                            strand = strand_pairing[strand]
                            nucs = utils.rev_seq(nucs)

                    except Exception as inst:
                        logger.debug("{0}: {1}".format(type(inst), inst.args))
                    
                    if 'N' not in nucs:
                        print >> pos_to_context, '\t'.join([gene, chrom + ':' + str(pos), strand, nucs, str(expr_dict[gene]), str(rep_dict[gene])])

    hin.close()
    pos_to_context.close()

    logger.info('finished getting trinucleotide context')

    # sort
    hout = open('merged_pos_to_context_sorted.txt', 'w')
    s_ret = subprocess.call(["sort", "-k2", "merged_pos_to_context.txt"], stdout = hout)
    hout.close()

    # sometimes, multiple genes share same genome as coding regions.
    #chr1    367658  368597  OR4F16(NM_001005277)    coding  +
    #chr1    367658  368597  OR4F29(NM_001005221)    coding  +
    # in these cases, we use higher expression.
    
    #concern: is there any genomic position which diffrent genes on different strands share?  --> we checked the posssibility and there was no position like that.

    rep_pos_df = pd.read_table("merged_pos_to_context_sorted.txt", sep='\t', names=('gene', 'pos', 'strand', 'trinucleotide', 'expression', 'replication_time'))
    del rep_pos_df['gene']
    rep_pos_df = rep_pos_df.drop_duplicates()
    rep_pos_df.to_csv("merged_pos_to_context_duplicate_dropped.txt", sep='\t', header =False, index=False)

    # top expression
    tmp_pos = ""
    tmp_exp = 0
    hout = open("merged_pos_to_context_duplicate_dropped2.txt", "w")
    with open("merged_pos_to_context_duplicate_dropped.txt", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if tmp_pos != F[0] and tmp_pos != "":
                print >> hout, '\t'.join([tmp_pos, tmp_strand, tmp_trinuc, tmp_exp, tmp_rep_time])
                tmp_pos = F[0]
                tmp_strand = F[1]
                tmp_trinuc = F[2]
                tmp_exp = F[3]
                tmp_rep_time = F[4]
                continue

            tmp_pos = F[0]
            tmp_strand = F[1]
            tmp_trinuc = F[2]
            if F[3] > tmp_exp:
                tmp_exp = F[3]
            tmp_rep_time = F[4]

        print >> hout, '\t'.join([tmp_pos, tmp_strand, tmp_trinuc, tmp_exp, tmp_rep_time])
    hin.close()
    hout.close()

    rep_pos_df = pd.read_table("merged_pos_to_context_duplicate_dropped2.txt", sep='\t', names=('pos', 'strand', 'trinucleotide', 'expression', 'replication_time'))
    rep_pos_df['replication_time'] = rep_pos_df['replication_time'].fillna(all_median_rep_time)
    rep_pos_df['exp_group'], exp_bins = pd.qcut(rep_pos_df['expression'].astype(float), 15, labels=['exp_class1', 'exp_class2', 'exp_class3', 'exp_class4', 'exp_class5', 
                                                                                                    'exp_class6', 'exp_class7', 'exp_class8', 'exp_class9', 'exp_class10',
                                                                                                    'exp_class11', 'exp_class12', 'exp_class13', 'exp_class14', 'exp_class15'
                                                                                                    ], retbins = True)
    rep_pos_df['rep_group'], rep_bins = pd.qcut(rep_pos_df['replication_time'].astype(float), 15, labels=['rep_class1', 'rep_class2', 'rep_class3', 'rep_class4', 'rep_class5', 
                                                                                                          'rep_class6', 'rep_class7', 'rep_class8', 'rep_class9', 'rep_class10',
                                                                                                          'rep_class11', 'rep_class12', 'rep_class13', 'rep_class14', 'rep_class15'
                                                                                                          ], retbins = True)
    rep_pos_df['class'] = rep_pos_df['exp_group'].str.cat(rep_pos_df['rep_group'], sep='_')

    rep_pos_df.to_csv("merged_pos_to_context_class.txt", sep = '\t', header =False, index=False)
    print "exp bins"
    print exp_bins
    print "rep bins"
    print rep_bins
    
    ######
    # output
    ######

    # output
    f = open('CCLE_avg_exp_output.txt', 'w')
    for k, v in sorted(expr_dict.items()):
        f.write('{0}\t{1}\n'.format(k, v))
    f.close()

    f2 = open('replication_time_output.txt', 'w')
    for k, v in sorted(rep_dict.items()):
        if v == 'NaN':
            rep_dict[k] = all_avg_rep_time
        f2.write('{0}\t{1}\n'.format(k, v))
    f2.close()


if __name__ == "__main__":
    db_preparation_main()