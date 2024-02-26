from typing import List
import random
from tqdm import tqdm
from motif_search.base_search import BaseMotifSearch
from utils.motif_utils import find_profile, find_kmer_prob, find_consensus
from utils.decorators import error_handler


class GibbsMotifSearch(BaseMotifSearch):
    """Find best Motif in DNA sequences based on gibbs sampling."""

    def __init__(
        self, 
        genes: List[str], 
        k: int,
        metric: str, 
        n_iter: int
    ):
        """Initialize Gibbs Motif Search algorithm.

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
        """Run Gibbs Motif Search algorithm.

        Returns:
            dict: respose with found motifs for each gene, their scores and consensus motif
        """    
        # set best motifs as random k-mers from each string in genes
        motifs = self.choose_random_motifs()
        best_motifs = motifs.copy()
        best_score = self.scoring_function(best_motifs)

        for _ in tqdm(range(self.n_iter)):
            i = random.randint(0, self.n_genes-1)                                 # choose random gene
            cut_motifs = [motif for idx, motif in enumerate(motifs) if idx != i]  # cut motifs without chosen gene
            profile = find_profile(cut_motifs)                                    # find new profile based on cut motifs
            updated_motif = self.sample_kmer(self.genes[i], profile)              # generate new motif based on updates profile
            motifs[i] = updated_motif                                             # update motifs with new motif
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
    
    def sample_kmer(self, gene: str, profile: dict) -> str:
        """Return a randomly generated k-mer from gene whose probabilities are generated from profile

        Args:
            gene (str): gene to sample k-mer from
            profile (dict): profile matrics based on cut motifs

        Returns:
            str: sampled k-mer
        """ 
        kmer_probs = dict()
        for i in range(0, self.last_index):
            kmer_probs[gene[i:i+self.k]] = find_kmer_prob(gene[i:i+self.k], profile)
        kmer_probs = self.normalize(kmer_probs)
        kmer = self.choose_kmer_with_die(kmer_probs)
        return kmer

    @staticmethod
    def normalize(probs: dict) -> dict:
        """Normalize probabilities: sum(probabilities) = 1.

        Args:
            probs (dict): input probabilities

        Returns:
            dict: normalized probabilities
        """
        probs_sum = sum(probs.values())
        normalized_probs = {key: value / probs_sum for key, value in probs.items()}
        return normalized_probs

    @staticmethod
    def choose_kmer_with_die(kmer_probs: dict) -> str:
        """Choose k-mer from dictionary based on probabilities of these k-mers by throwing a die

        Args:
            kmer_probs (dict): probabilities of each k-mer in gene

        Returns:
            str: choosen k-mer
        """
        random_value = random.uniform(0, 1)
        total_prob = 0
        for kmer, prob in kmer_probs.items():
            total_prob += prob
            if random_value < total_prob:
                return kmer