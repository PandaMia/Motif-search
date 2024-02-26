from typing import List
from motif_search.base_search import BaseMotifSearch
from utils.motif_utils import find_profile, find_most_probable_kmer, find_consensus
from utils.decorators import error_handler


class GreedyMotifSearch(BaseMotifSearch):
    """Greedy find best Motif in DNA sequences."""

    def __init__(
        self, 
        genes: List[str], 
        k: int,
        metric: str 
    ):
        """Initialize Greedy Motif Search algorithm.

        Args:
            genes (List[str]): list of genes: genes[str]
            k (int): length of the motif
            metric (str): metric to evaluate found motifs
        """
        super().__init__(genes, k, metric)  # initialize BaseMotifSearch

    @error_handler
    def run_search(self) -> dict:
        """Run Greedy Motif Search algorithm.

        Returns:
            dict: respose with found motifs for each gene, their scores and consensus motif
        """
        # set best_motifs as first k-mer from each string in genes
        best_motifs = []
        for i in range(0, self.n_genes):
            best_motifs.append(self.genes[i][0:self.k])

        # iterate over each k-mer from the first row of genes
        for i in range(self.gene_len - self.k + 1):
            motifs = []
            motifs.append(self.genes[0][i:i+self.k])   # set current k-mer as first row in Motifs
            for j in range(1, self.n_genes):           # iterate over other rows of genes
                profile = find_profile(motifs[0:j])    # set current profile based on chosen rows for Motif (on the first step this is only one row)
                motifs.append(find_most_probable_kmer(self.genes[j], self.k, profile))  # add most probable kmer from the next row of genes to Motifs 

            # after building whole Motif matrix compare it with previous Best Motif
            if self.scoring_function(motifs) < self.scoring_function(best_motifs):
                best_motifs = motifs.copy()
        
        # get scores for the best motifs and consensus motif
        scores = self.evaluate_best_motifs(best_motifs)
        consensus = find_consensus(best_motifs)

        response = {
            "best_motifs": best_motifs,
            "scores": scores,
            "consensus": consensus
        }
        
        return response