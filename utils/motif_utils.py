from typing import List


PSEUDOCOUNT: int = 1


def count_nucleotides(motifs: list) -> dict:
    """Count the number of nucleotides on each position in motifs with pseudocounts.

    Args:
        motifs (list): motifs matrix

    Returns:
        dict: dictionary with counts of nucleotides: counter[nucleotide][position]
    """
    t = len(motifs)
    k = len(motifs[0])

    count = {}
    for nucleotide in "ACGT":
        count[nucleotide] = [PSEUDOCOUNT] * k

    for row in range(t):
        for position in range(k):
            nucleotide = motifs[row][position]
            count[nucleotide][position] += 1
    return count


def find_consensus(motifs: list) -> str:
    """Find most popular nucleotide in each column of the motif matrix.

    Args:
        motifs (list): motifs matrix

    Returns:
        str: consensus motif
    """
    # 
    k = len(motifs[0])
    nucleotide_counter = count_nucleotides(motifs)
    consensus = ""
    for j in range(k):
        m = 0
        most_frequent_nucleotide = ""
        for nucleotide in "ACGT":
            if nucleotide_counter[nucleotide][j] > m:
                m = nucleotide_counter[nucleotide][j]
                most_frequent_nucleotide = nucleotide
        consensus += most_frequent_nucleotide
    return consensus


def find_profile(motifs: list) -> dict:
    """Find frequency of nucleotides on each position in motifs with pseudocounts

    Args:
        motifs (list): motifs matrix

    Returns:
        dict: dictionary with frequencies of nucleotides: profile[nucleotide][position]
    """
    # 't' is the number of motifs + number of psuedocounts
    t = len(motifs) + len("ACGT") * PSEUDOCOUNT
    k = len(motifs[0])

    profile = {}
    for nucleotide in "ACGT":
        profile[nucleotide] = [0] * k
    
    count = count_nucleotides(motifs)
    for nucleotide in count:
        for position in range(len(count[nucleotide])):
            profile[nucleotide][position] = count[nucleotide][position] / t
    return profile


def find_kmer_prob(kmer: str, profile: dict) -> float:
    """Calculate probability of appearing specific k-mer with given Profile matrix

    Args:
        kmer (str): given k-mer
        profile (dict): dictionary with probabilities to see each nucleotide on each position: profile[nucleotide][position]

    Returns:
        float: probability of appearing k-mer
    """
    prob = 1
    for i in range(len(kmer)):
        prob *= profile[kmer[i]][i]
    return prob


def find_most_probable_kmer(gene: str, k: int, profile: dict) -> str:
    """Find the most probable k-mer in the gene with given Profile matrix

    Args:
        gene (str): string to search for k-mer
        k (int): length of k-mer
        profile (dict): dictionary with probabilities to see each nucleotide on each position: profile[nucleotide][position]

    Returns:
        str: the most probable k-mer
    """
    # 
    max_prob = 0
    most_prob_kmer = gene[:k]
    for i in range(len(gene) - k + 1):
        kmer = gene[i:i+k]
        prob = find_kmer_prob(kmer, profile)
        if prob > max_prob:
            most_prob_kmer = kmer
            max_prob = prob
    return most_prob_kmer


def choose_motifs(genes: List[str], k: int, profile: dict) -> List[str]:
    """Find the most probable k-mer for each gene with given Profile matrix

    Args:
        genes (List[str]): list of genes: genes[str]
        k (int): length of k-mer
        profile (dict): dictionary with probabilities to see each nucleotide on each position: profile[nucleotide][position]

    Returns:
        List[str]: list of most probable k-mers
    """
    motifs = [find_most_probable_kmer(gene, k, profile) for gene in genes]
    return motifs