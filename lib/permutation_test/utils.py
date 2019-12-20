#!/usr/bin/env python2

import os
import subprocess
import argparse
import logging
import warnings
import sys
from functools import wraps
import datetime
import pkg_resources
import gzip

logger = logging.getLogger(__name__)  # module logger

# reference: probabilistic2020 (https://github.com/KarchinLab/probabilistic2020  under Apache License 2.0)

# global dictionary specifying base pairing
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


def check_valid_nuc(nuc):
    valid_nucs = ['A', 'C', 'T', 'G', 'N']
    is_valid = nuc in valid_nucs
    return is_valid


def reverse_complement(seq):
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base_pairing[s] for s in rev_seq])
    return rev_comp_seq


def start_logging(log_file='stdout', log_level='INFO', verbose=False):

    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO

    # ignore warnings if not in debug
    if log_level.upper() != 'DEBUG':
        warnings.filterwarnings('ignore')

    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

    # logging to stdout
    root = logging.getLogger()
    root.setLevel(lvl)
    stdout_stream = logging.StreamHandler(sys.stdout)
    stdout_stream.setLevel(lvl)
    formatter = logging.Formatter(myformat)
    stdout_stream.setFormatter(formatter)
    root.addHandler(stdout_stream)
    root.propagate = True
