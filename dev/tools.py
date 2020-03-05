#!/usr/bin/env python3

#############################################################
#
# tools.py
#
# Description : set of tools to perform sequence checking
#
# Authors     : Samuel Miravet-Verde, Anas Gharrab, Marc Weber,
#               Maria Lluch-Senar, Luis Serrano
#
#############################################################

import json
import os

import pandas as pd
import regex as re
from Bio import SeqIO
from Bio.Alphabet import IUPAC, _verify_alphabet, generic_dna
from Bio.Seq import Seq
from RNA import fold
from Bio.SeqUtils import GC

from ribosome_binding_site import RBS_CANONICAL
from ribosome_binding_site import compute_RBS_affinity


###
# I/O functions
###

def is_dna_seq_valid(value, op=True):
    """ Validates DNA sequence. If op is True check IUPACUnambiguousDNA library if False IUPACAmbiguousDNA"""
    return _verify_alphabet(Seq(value, IUPAC.IUPACUnambiguousDNA())) if op else _verify_alphabet(Seq(value, IUPAC.IUPACAmbiguousDNA()))

def seq_translator(seq, to='DNA'):
    """ DNA/RNA to Protein sequence """
    generic = generic_dna if to == 'DNA' else generic_rna
    return Seq(seq, generic).translate()

def load_sequence(sequence):
    """ If sequence is a fasta file, return sequence """
    if os.path.isfile(sequence):
        handle = open(sequence, 'rU')
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
        handle.close()
    return sequence


def check_sequence_iterator(sequence, n, circular=False):
    """
    Given a <sequence>, an integer <n> corresponding to the desired window size and a boolean for <circular>
    Returns [sequence, limit] adjusted for circular construct or linear
    """
    l = len(sequence)
    if l < n:
        raise ValueError('Sequence to evaluate shorter than evaluation matrix')
    else:
        if circular:
            # Require to check plasmid that are circular
            sequence = sequence + sequence
            limit = l
        else:
            limit = l - n + 1
        return [sequence, limit]


def check_outhandle(d, outhandle='json', window_len=0):
    """ Converts <d> in shape {<position>:<score>} to [{"start":"<position>", "score":"<score>"}] """
    if outhandle == 'json':
        return json.dumps([{'start':k, 'end':k+window_len, 'raw_score': round(v, 2), 'norm_score': round(v, 2)} for k, v in d.items()])
    else:
        return d


def load_matrix(inFile, residue_type='DNA', k=1, indexed=True):
    """
    Loads a matrix from text file, if nucl, matrix is assumed to have the shape mxn
    where:
        m is the number of rows, i.e. nucleotide bases or amino acids
        n is the number of columns, i.e. positions of the motif to be evaluated
      i  i+1  i+2, ..., i+n
    A
    C
    G
    T
    indexes can be or not included, if not included, alphabetical order is assumed.
    In the case of k-mers, parameter k can be set to > 1. Example: DNA triplets (codons).
    Returns: the matrix in compatible format
    """
    if os.path.isfile(inFile):
        if indexed:
            matrix = pd.read_csv(inFile, sep='\t', index_col=0)  # If matrix format changes, update it here
        else:
            matrix = pd.read_csv(inFile, sep='\t', header=None)  # If matrix format changes, update it here
        # Check if properly formated
        if ((residue_type in ['DNA', 'RNA'] and matrix.shape[0] != 4 ** k) or
                (residue_type not in ['DNA', 'RNA'] and matrix.shape[0] != 20 ** k)):
            raise ValueError('Matrix has {} rows, expected {} for {}.'.format(matrix.shape[0], 4 ** k, residue_type))
        else:
            return matrix
    else:
        raise ValueError('Matrix not found, please provide a compatible input path')

###
# IUPAC PATTERN MATCHING
###

def seq2regex(seq):
    IUPAC = {"A" : "A",
             "C" : "C",
             "G" : "G",
             "T" : "T",
             "U" : "U",
             "R" : "[GA]",
             "Y" : "[TC]",
             "K" : "[GT]",
             "M" : "[AC]",
             "S" : "[GC]",
             "W" : "[AT]",
             "B" : "[CGT]",
             "D" : "[AGT]",
             "H" : "[ACT]",
             "V" : "[ACG]",
             "N" : "."}
    return ''.join([IUPAC.get(i) for i in seq])

def match_sequence(subsequence, sequence):
    """
    Given an ambiguous <subsequence> and a DNA <sequence>,
    Returns an array with the positions where <subsequence>
    matches <sequence>.
    Example:
        match_sequence('ABMHA', 'NNNNAGCTANNNN') --> [4]
        ABMHA should match A[C|G|T][A|C][A|C|T]A
    """
    return [(m.start(0)+1, m.end(0)) for m in re.finditer(seq2regex(subsequence), sequence, overlapped=True)]

###
# Specific evaluation functions
###
def matrix_scoring(sequence, matrix, circular=False, residue_type='DNA',
                   indexed=True, outhandle='json', standardize=True):
    """
    Scoring function
    """
    # Check matrix and load it if required
    print('\n\n--> Matrix scoring started for {}\n\tMatrix scoring: Loading scoring matrix...'.format(matrix))
    if type(matrix) == str:
        matrix = load_matrix(matrix, residue_type=residue_type)
    # Iterate by windows and score
    n = matrix.shape[1]
    sequence, limit = check_sequence_iterator(sequence, n, circular)
    print('\tMatrix scoring: Scoring sequence...')
    rs = {i + 1: sum([matrix.at[sequence[i:i + n][subindex], str(subindex + 1)] for subindex in range(n)]) for i in
          range(0, limit)}
    print('\tMatrix scoring: Returning results...\n--> Matrix scoring finished.\n\n')
    return check_outhandle(rs, outhandle, n)


def RNAstructure_scoring(sequence, n=20, circular=False, residue_type='DNA', outhandle='json', standardize=True):
    """
    Given a <sequence>, an integer <n> corresponding to the desired window size,
    Returns a <outhandle> object with the secondary structure energy for the different windows
    """
    # Check matrix and load it if required
    print('\n\n--> RNA structure scoring started for window size={}'.format(n))
    # Iterate by windows and score
    sequence, limit = check_sequence_iterator(sequence, n, circular)
    print('\tRNA structure scoring: Scoring sequence...')
    rs = {i + 1: fold(sequence[i:i + n])[1] for i in range(0, limit)}
    print('\tRNA structure scoring: Returning results...\n--> RNA structure scoring finished.\n\n')
    return check_outhandle(rs, outhandle, n)

def GC_scoring(sequence, n=20, circular=False, residue_type='DNA', outhandle='json', standardize=True):
    """
    Given a <sequence>, an integer <n> corresponding to the desired window size,
    Returns a <outhandle> object with the GC content for the different windows
    """
    # Check matrix and load it if required
    print('\n\n--> scoring started for window size={}'.format(n))
    # Iterate by windows and score
    sequence, limit = check_sequence_iterator(sequence, n, circular)
    print('\tGC scoring: Scoring sequence...')
    rs = {i + 1: GC(sequence[i:i + n]) for i in range(0, limit)}
    print('\tGC scoring: Returning results...\n--> GC scoring finished.\n\n')
    return check_outhandle(rs, outhandle, n)


def extract_codons_list(seq, frame=0, checkLengthMultipleOf3=False, frameFromEnd=False):
    if len(seq) % 3 != 0 and checkLengthMultipleOf3:
        print("ERROR: seq length is not multiple of 3.")
        return None

    l = len(seq)
    if frameFromEnd:
        frame1 = ((l % 3) + frame) % 3
    else:
        frame1 = frame
    codonList = (seq[3 * n + frame1: 3 * n + frame1 + 3] for n in range(0, int((l - frame1) / 3)))
    return codonList


def codon_adaptation_scoring(sequence, matrix, circular=False, residue_type='DNA', n=1,
                             indexed=True, outhandle='json', standardize=True, verbose=1):
    """
    Returns a score reflecting local codon adaptation, computed based on codon usage, averaged
    over a window of size n.
    In this case the matrix should be a 1-column table with codon index and codon score values.
    We use the w_ij coefficient of the normalized relative codon usage, such that the most frequent
    codon in each synonymous group has a score of 1.
    """
    # Check matrix and load it if required
    print('\n\n--> Matrix scoring started for {}\n\tMatrix scoring: Loading scoring matrix...'.format(matrix))
    if type(matrix) == str:
        matrix = load_matrix(matrix, residue_type=residue_type, k=3)
        if verbose >= 2:
            print("matrix\n", matrix)

    sequence, limit = check_sequence_iterator(sequence, n, circular)
    if verbose >= 1: print('\tCodon adaptation scoring: Scoring sequence...')
    codons = extract_codons_list(sequence, checkLengthMultipleOf3=False)
    # Averaged rolling window
    rs = pd.Series(codons).map(matrix.iloc[:, 0].to_dict()).rolling(window=n).mean()
    # convert codon position to nucleotide position, repeat codon score
    # for each nucleotide position in the codon
    rs0 = rs.copy()
    rs0.index = 3 * rs0.index
    rs1 = rs.copy()
    rs1.index = 3 * rs1.index + 1
    rs2 = rs.copy()
    rs2.index = 3 * rs2.index + 2
    rs = pd.concat([rs0, rs1, rs2]).sort_index()
    # Convert to 1-based index
    rs.index = rs.index + 1
    rs = rs.to_dict()
    if verbose >= 1: print('\tCodon adaptation scoring: Returning results...\n--> Matrix scoring finished.\n\n')
    return check_outhandle(rs, outhandle, n)


def RBS_scoring(sequence, motif=RBS_CANONICAL, circular=False, residue_type='DNA', start_codons=['ATG', 'GTG', 'TTG'],
                indexed=True, outhandle='json', standardize=True):

    """Detects the presence of ribosome binding site motif upstream of any stop codons in the sequence.
    Returns a score for each RBS corresponding to its hybridization energy to the anti-motif.

    Note: The rs dictionary only contains an entry for the positions with a detected RBS.
    """

    pattern = '(' + '|'.join(start_codons) + ')'
    rs = dict()
    RBSs = []
    for m in re.finditer(pattern, sequence):
        start = m.start()
        energy, rbs_start, fold = \
            compute_RBS_affinity(motif, sequence, start_codon_position=start, spacer_range=[5, 11],
                                 include_internal_motif_around_start_codon=False, verbose=0)
        rbs_start = int(rbs_start)
        # very small threshold to consider the motif as a RBS
        if energy < -0.5:
            for i in range(len(motif)):
                rs[rbs_start + i] = energy
            RBSs.append((rbs_start, energy))

    print('\RBS_scoring: finished.\n\n')
    return check_outhandle(rs, outhandle, len(motif))   # window len not defined


def fixed_matrix_scoring(sequence, matrix, circular=False, residue_type='DNA', fixed_sequences = ['ATG', 'GTG', 'TTG'], mode=1,
                         indexed=True, outhandle='json', standardize=True):
    """
    matrix scoring function, that runs only if any of the <fixed_sequences> is present in the end of each subsequence in <sequence> (mode 1)
    or at the beginning (mode!=1). The fixed sequences are not considered in the evaluation.

    fixed_sequences is expected to have sequences all with the same length
    """
    # Check matrix and load it if required
    print('\n\n--> Matrix scoring started for {}\n\tMatrix scoring: Loading scoring matrix...'.format(matrix))
    if type(matrix) == str:
        matrix = load_matrix(matrix, residue_type=residue_type)

    # Iterate by windows and score
    extra = max([len(i) for i in fixed_sequences])
    n = matrix.shape[1] # THIS IS DIFFERENT
    sequence, limit = check_sequence_iterator(sequence, n+extra, circular)
    print('\tMatrix scoring: Scoring sequence...')

    if len(fixed_sequences)==0:
        raise ValueError('fixed_sequences empty')
    else:
        if type(fixed_sequences)==str:
            fixed_sequences = [fixed_sequences]

    rs = {}
    for i in range(0, limit):
        subseq = sequence[i:i+n+extra]
        if mode==1:
            if any([fixseq==subseq[-extra:] for fixseq in fixed_sequences]):
                # print(i, subseq)
                subseq = subseq[:-extra]
                # print(i, subseq)
                rs[i+1]=sum([matrix.at[subseq[subindex], str(subindex + 1)] for subindex in range(n)])
            else:
                rs[i+1]=0
        else:
            if any([fixseq==subseq[:extra] for fixseq in fixed_sequences]):
                # print(i, subseq)
                subseq = subseq[extra:]
                # print(i, subseq)
                rs[i+1]=sum([matrix.at[subseq[subindex], str(subindex + 1)] for subindex in range(n)])
            else:
                rs[i+1]=0

    print('\tMatrix scoring: Returning results...\n--> Matrix scoring finished.\n\n')
    return check_outhandle(rs, outhandle, n)

###
# Main checker
###


def checker(sequence,
            elements='all', matrix_dict=None,
            codon_table=4,
            circular=True, residue_type='DNA',
            indexed=True, outhandle='json', standardize=True, verbose=0):
    # Parse and check elements to explore
    default_elements = set(['promoter', 'terminator', 'utr5', 'codon_adaptation', 'NTPi', 'iRNA',
                            'restriction_sites', 'RBS', 'toxic', 'alternative_start', 'RNA_structure20',
                            'RNA_structure60', 'GC20'])
    if elements == 'all':
        elements = default_elements
    elif type(elements) == str:
        elements = set([elements])
    elif type(elements) != set:
        elements = set(elements)
    # Take only considered elements
    if len(elements.difference(default_elements)) != 0:
        print('{} checker not implemented, it will not be calculated.'.format(
            list(elements.difference(default_elements))[0]))
        elements = elements.intersection(default_elements)

    # Generate all results
    rs = {}
    for element in elements:
        if element in matrix_dict:
            if element == 'codon_adaptation':
                rs[element] = codon_adaptation_scoring(sequence, matrix_dict[element], n=1,
                                                       circular=circular, residue_type=residue_type,
                                                       indexed=indexed, outhandle=outhandle, standardize=standardize,
                                                       verbose=verbose)
            elif element == 'utr5':
                rs[element] = fixed_matrix_scoring(sequence, matrix_dict[element], circular=circular,
                                                   residue_type='DNA', fixed_sequences = ['ATG', 'GTG', 'TTG'], mode=1,
                                                   indexed=indexed, outhandle=outhandle, standardize=standardize)
            else:
                rs[element] = matrix_scoring(sequence, matrix_dict[element], circular=circular,
                                             residue_type=residue_type,
                                             indexed=indexed, outhandle=outhandle, standardize=standardize)
        elif element == 'RNA_structure20':
            rs[element] = RNAstructure_scoring(sequence, n=20, circular=circular, residue_type=residue_type,
                                               outhandle=outhandle, standardize=standardize)
        elif element == 'RNA_structure60':
            rs[element] = RNAstructure_scoring(sequence, n=60, circular=circular, residue_type=residue_type,
                                               outhandle=outhandle, standardize=standardize)
        elif element == 'GC20':
            rs[element] = GC_scoring(sequence, n=20, circular=circular, residue_type=residue_type,
                                     outhandle=outhandle, standardize=standardize)
        # elif element == 'RBS':
        #     rs[element] = RBS_scoring(sequence, circular=circular, outhandle=outhandle)
    return rs   # TODO second dictionary is expected to be the normalized dictionary
