#!/usr/bin/env python3

#############################################################
#
# bioroboost.py
#
# Author      : Miravet-Verde, Samuel
# Co-authors  : Anas Gharrab, Marc Weber
# Supervision : Maria Lluch-Senar, Luis Serrano
#
#############################################################

import sys, os
import argparse
import tools as brt

###
# GLOBALS
###

# residues = {'DNA':['A','C','G','T'], 'RNA':['A','C','G','U'], 'PRO':['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']}

###
# MAIN FUNCTION
###

def main(mode=1, **kwargs):
    """
    modes:
        1-scoring
        2-training
        3-consensus
    """
    if mode=='1' or mode==1 or mode=='scoring':
        if all(i in kwargs for i in ['sequence', 'matrix']):
            result = brt.matrix_scoring(sequence=brt.load_sequence(kwargs['sequence']), matrix=kwargs['matrix'])
            print(result)
            return result
        else:
            sys.exit('No sequence or matrix found for scoring function')

    elif mode==2 or mode=='training':
        pass
    else:
        pass

###
# ARGPARSE
###

if '__main__'==__name__:

    # CMD PARSER
    parser = argparse.ArgumentParser(description = "bioroboost allows the detection of subsequences in a sequence that could prevent or complicate the application of that sequence in synthetic biology tasks.")
    parser.add_argument('-f', '--function',
                        dest="mode",
                        action="store",
                        default='scoring',
                        required=True,
                        help="Mode of function:\n\t- 1. scoring\n\t- 2. training\n\t- 3. consensus")

    parser.add_argument('-s', '--sequence',
                        dest="sequence",
                        action="store",
                        required=True,
                        type=str,
                        help="User sequence")

    parser.add_argument('-m', '--matrix',
                        dest="matrix",
                        action="store",
                        required=True,
                        type=str,
                        help="Element matrix to use")
    options = parser.parse_args()

    # Run
    main(mode=options.mode, sequence=options.sequence, matrix=options.matrix)
