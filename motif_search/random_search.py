from typing import List
import random
from tqdm import tqdm
from motif_search.base_search import BaseMotifSearch
from utils.motif_utils import find_profile, choose_motifs, find_consensus
from utils.decorators import error_handler


class RandomMotifSearch(BaseMotifSearch):
    """Random find best Motif in DNA sequences."""

    def __init__(
        self, 
        genes: List[str], 
        k: int,
        metric: str, 
        n_iter: int
    ):
        """Initialize Random Motif Search algorithm.

        Args:
            genes (List[str]): list of genes: genes[str]
            k (int): length of the motif
            metric (str): metric to evaluate found motifs
            n_iter (int): number of iterations to search motifs
        """
        super().__init__(genes, k, metric)        # initialize BaseMotifSearch
        self.n_iter = n_iter                      # number of epochs
        self.last_index = self.gene_len - self.k  # last possible index for motif start

    @error_handler
    def run_search(self) -> dict:
        """Run Random Motif Search algorithm.

        Returns:
            dict: respose with found motifs for each gene, their scores and consensus motif
        """
        best_motifs = None
        best_score = float("inf")

        # run randomized search N times
        for _ in tqdm(range(self.n_iter)):
            motifs = self.run_epoch()
            score = self.scoring_function(motifs)
            if score < best_score:
                best_score = score
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

    def run_epoch(self) -> List[str]:
        """Run epoch of random search algorithm.

        Returns:
            List[str]: best motifs for the epoch
        """
        # set best motifs as random k-mers from each string in genes
        motifs = self.choose_random_motifs()
        best_motifs = motifs.copy()

        best_score = self.scoring_function(best_motifs)
        while True:
            profile = find_profile(motifs)
            motifs = choose_motifs(self.genes, self.k, profile)
            score = self.scoring_function(motifs)
            if score < best_score:
                best_score = score
                best_motifs = motifs.copy()
            else:
                return best_motifs
    
    def choose_random_motifs(self) -> List[str]:
        """Choose random motifs from each gene.

        Returns:
            List[str]: random motifs
        """
        random_motifs = []
        for i in range(self.n_genes):
            index = random.randint(0, self.last_index)
            motif = self.genes[i][index: index + self.k]
            random_motifs.append(motif)
        return random_motifs