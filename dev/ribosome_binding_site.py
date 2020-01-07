#!/usr/bin/env python3

import numpy as np
from Bio.Seq import Seq
import RNA


RBS_CANONICAL = 'AGGAGG'


def compute_best_motif_hybridization_energy_in_seq(motif, seq, reverse_complement_motif=True, scan=True, verbose=1):
    """
    Computes the hybridization energy of the motif (e.g. Shine-Dalgarno AGGAGG) for every subsequence
    in `seq` and returns the energy and position of the lowest energy.

    If `reverse_complement_motif` is False, the `motif` argument should be already the complementary motif
    (e.g. Shine-Dalgarno CCTCCT).
    """
    if type(motif) is not Seq:
        motif_bio = Seq(motif)
    else:
        motif_bio = motif
    if reverse_complement_motif:
        motif_seq = str(motif_bio.reverse_complement())
    else:
        motif_seq = str(motif_bio)
    if scan:
        if verbose >= 2: print("motif_seq", motif_seq)
        l = len(motif_seq)
        if len(seq) < l:
            subseq_list = [seq]
        else:
            subseq_list = [seq[i:i+l] for i in range(len(seq) - l + 1)]
        if verbose >= 2: print("subseq_list", subseq_list)
        energy_fold_list = [RNA.cofold(motif_seq + '&' + subseq) for subseq in subseq_list]
        energy_list = [e[1] for e in energy_fold_list]
        if verbose >= 2: print("energy_fold_list", energy_fold_list)
        energy_index = np.argmin(energy_list)
        energy_tuple = (energy_list[energy_index], energy_index, energy_fold_list[energy_index][0])
    else:
        energy_tuple = (RNA.cofold(motif_seq + '&' + seq)[1], 0)
    return energy_tuple


def compute_RBS_affinity(motif, seq, start_codon_position, UTR5p_position=0, spacer_range=[5, 11],
                         include_internal_motif_around_start_codon=False, verbose=1):
    """Computes the strongest RBS upstream of the start codon, within a spacer distance range.
    
    Example:
    ```
sequence = ('CCACTCGA'    # not transcribed
            'AGACGGGACCGCCGGCGAGCAGGTGCACAAAGCC'    # 5'UTR
            'ATGAAGCGCTACGCCCTGGTGCCCGGCACCATCGCCTTTACCGACGCACATATCGAGGTGGACATTACCTACGCCGAGTACTTCTAA')    # CDS
compute_RBS_affinity('AGGAGG', sequence, start_codon_position=42, UTR5p_position=8, spacer_range=[5, 11],
                     include_internal_motif_around_start_codon=False, verbose=2)
    ```
    """
    assert spacer_range[0] >= 0
    assert spacer_range[1] >= 0
    region = ''
    region1 = max(-(len(motif) + spacer_range[1]) + start_codon_position, UTR5p_position)
    region2 = -(spacer_range[0] + 1) + start_codon_position
    if include_internal_motif_around_start_codon:
        # Examine region around start codon, up to 6 nt downstream of the ATG
        region2 = start_codon_position + 9
    if region1 < region2:
        region += seq[region1: region2 + 1]
        if verbose >= 2:
            print("###")
            print("motif", motif, "region1", region1, "region2", region2)
            print("spacer_range", spacer_range)
            if region2 <= start_codon_position:
                print('|'.join([seq[:UTR5p_position],
                                seq[UTR5p_position:region1], seq[region1:region2 + 1],
                                seq[region2 + 1:start_codon_position],
                                seq[start_codon_position:start_codon_position+3]]))
            else:
                print('|'.join([seq[:UTR5p_position],
                                seq[:UTR5p_position:region1], seq[region1:start_codon_position],
                                seq[start_codon_position:region2 + 1]]))
        energy, energy_index, fold = compute_best_motif_hybridization_energy_in_seq(motif, region, verbose=verbose)
        energy_index = energy_index + region1
        return energy, energy_index, fold
    else:
        return None, None, None
