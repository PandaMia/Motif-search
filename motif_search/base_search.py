from typing import List
from config.config import SCORING_FUNCTIONS
from utils.evaluation import calc_entropy_score, calc_hamming_distance
from utils.motif_utils import find_consensus
from utils.decorators import error_handler


class BaseMotifSearch():
    """Base class for motif search algorithms."""

    def __init__(
        self, 
        genes: List[str], 
        k: int,
        metric: str
    ):
        """Initialize Base motif search class.

        Args:
            genes (List[str]): list of genes: genes[str]
            k (int): length of the motif
            metric (str): metric to evaluate found motifs
        """
        self.genes = genes                                 # list of genes to search motifs
        self.k = k                                         # length of the motif
        self.n_genes = len(genes)                          # number of genes
        self.gene_len = len(genes[0])                      # length of the single gene
        self.scoring_function = SCORING_FUNCTIONS[metric]  # function to evaluate found motifs

    @error_handler
    def run_search(self) -> dict:
        """Run motif search algorithm.

        Returns:
            dict: respose with found motifs for each gene, their scores and consensus motif
        """
        # set best motifs as first k-mer from each string in genes
        best_motifs = []
        for i in range(0, self.n_genes):
            best_motifs.append(self.genes[i][0:self.k])
        
        # get scores for the best motifs and consensus motif
        scores = self.evaluate_best_motifs(best_motifs)
        consensus = find_consensus(best_motifs)

        response = {
            "best_motifs": best_motifs,
            "scores": scores,
            "consensus": consensus
        }
        return response
    
    def evaluate_best_motifs(self, best_motifs: List[str]) -> dict:
        """Evaluate found motifs.

        Args:
            best_motifs (List[str]): found best motifs

        Returns:
            dict: scores of the found motifs
        """
        scores = dict()
        scores["hamming_distance"] = calc_hamming_distance(best_motifs)
        scores["entropy_score"] = round(calc_entropy_score(best_motifs), 3)
        return scores