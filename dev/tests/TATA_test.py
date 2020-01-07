#!/usr/bin/env python3

#############################################################
#
# TATA_test.py
#
# Author : Miravet-Verde, Samuel
# Last update :  02-July-2019
#
#############################################################

import sys, os
from Bio import motifs, AlignIO
from Bio.Seq import Seq
#from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio.Alphabet import IUPAC, Gapped

#def run_clustalw(infile, outfile):
#    """
#    Run a basic clustal alignment and return a biopython alignment object
#    TO IMPLEMENT
#    Gap open penalty: -gapopen
#    Gap extension penalty: -gapext
#    no end gap(yes, no): -endgaps
#    gap distance: -gapdist
#    weight matrix(blosum, pam etc): -matrix ["BLOSUM", "PAM", "GONNET", "ID"]
#    type (DNA,protein): -type
#    """
#    cline = ClustalwCommandline("clustalw", infile=infile, outfile=outfile,
#                                seqnos="ON", gapopen=2, gapext=0.5)
#    os.system(str(cline))

def run_muscle(infile, outfile):
    """ Muscle alignment """
    cline = 'muscle -in '+infile+' -out '+outfile+' -clwstrict -quiet -log '+outfile+'.log'
    os.system(cline)

def generate_motif(infile, keep_intermediates=False):

    # Define outfile and alphabet
    outfile = os.path.splitext(infile)[0]+'.aln'
    # run alignment
    #run_clustalw(infile, outfile)
    run_muscle(infile, outfile)

    #Loading the msa and building a sequence motif
    alphabet = Gapped(IUPAC.unambiguous_dna)
    alignment_file = AlignIO.read(outfile,'clustal',
                              alphabet=alphabet)
    instances = []
    for sequence in alignment_file:
        instances.append(sequence.seq)
    m = motifs.create(instances, alphabet=alphabet)

    return m, alignment_file
    # print(m.pwm())

def background_probabilities(GC_content=0.4):
    """
    Given a GC content (1 base)
    Returns a dictionary with the background probability for each DNA nucleotide
    """
    return {'A':(1-GC_content)/2.0, 'C':GC_content/2.0, 'G':GC_content/2.0, 'T':(1-GC_content)/2.0, '-':0.0}

def evaluate_sequence(sequence, sequences_file='./promoters_TATAAT.fa', GC_content=0.4):
    """
    Given a string <sequence> and an <alignment_file> in clustal or fasta format
    Returns the pssm score of the <sequence>
    """
    t, aln= generate_motif(sequences_file)
    pwm = t.counts.normalize(pseudocounts=0.5)
    background = background_probabilities(GC_content)
    pssm = pwm.log_odds(background)
    pssm.alphabet =  IUPAC.unambiguous_dna
    _ = pssm.pop('-')
    return pssm.calculate(Seq(sequence, alphabet=IUPAC.unambiguous_dna))

if __name__ == '__main__':
    if len(sys.argv)==2:
        print(evaluate_sequence(sequence=sys.argv[1], sequences_file='./promoters_TATAAT.fa'))
    elif len(sys.argv)==3:
        print(evaluate_sequence(sequence=sys.argv[1], sequences_file=sys.argv[2]))
    else:
        print(evaluate_sequence(sequence='TAAAAA', sequences_file='./promoters_TATAAT.fa'))
        print(evaluate_sequence(sequence='TATAAT', sequences_file='./promoters_TATAAT.fa'))
