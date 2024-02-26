import numpy as np
from utils.motif_utils import find_consensus, count_nucleotides


def calc_hamming_distance(motifs: list) -> int:
    """Calculate score for Motifs matrix by comparing to consensus motif.

    Args:
        motifs (list): motifs matrix

    Returns:
        int: score
    """
    t = len(motifs)
    k = len(motifs[0])
    score = 0
    consensus = find_consensus(motifs)
    for i in range(t):
        for j in range(k):
            if motifs[i][j] != consensus[j]:
                score += 1
    return score


def calc_entropy_score(motifs: list) -> float:
    """Calculate entropy by nucleotides for given Motifs matrix.

    Args:
        motifs (list): motifs matrix

    Returns:
        float: entropy score
    """
    nucleotide_counter = count_nucleotides(motifs)
    matrix = np.array(list(nucleotide_counter.values()))
    probs = np.divide(matrix, np.sum(matrix, axis=0))
    score = -np.sum(probs * np.log2(
        probs, 
        out=np.zeros_like(probs), 
        where=(probs != 0)
    ))
    return score