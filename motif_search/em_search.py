from typing import List
import random

import numpy as np
from tqdm import tqdm

from motif_search.base_search import BaseMotifSearch
from utils.decorators import error_handler
from utils.motif_utils import PSEUDOCOUNT, find_consensus


class EMMotifSearch(BaseMotifSearch):
    """Expectation-maximization motif search with one motif per sequence."""

    def __init__(
        self,
        genes: List[str],
        k: int,
        metric: str,
        n_iter: int,
        n_starts: int,
    ):
        super().__init__(genes, k, metric)
        self.n_iter = max(1, n_iter or 1)
        self.n_starts = max(1, n_starts or 1)
        self.last_indexes = [len(gene) - self.k for gene in self.genes]
        self.convergence_tol = 1e-6

    @error_handler
    def run_search(self) -> dict:
        best_motifs = None
        best_score = float("inf")
        best_likelihood = float("-inf")

        for _ in tqdm(range(self.n_starts)):
            motifs, likelihood = self.run_restart()
            score = self.scoring_function(motifs)
            if score < best_score or (score == best_score and likelihood > best_likelihood):
                best_motifs = motifs.copy()
                best_score = score
                best_likelihood = likelihood

        scores = self.evaluate_best_motifs(best_motifs)
        consensus = find_consensus(best_motifs)

        response = {
            "best_motifs": best_motifs,
            "scores": scores,
            "consensus": consensus,
        }
        return response

    def run_restart(self):
        starts = self.choose_random_starts()
        motifs = self.extract_motifs(starts)
        profile = self.build_profile_from_motifs(motifs)
        previous_likelihood = float("-inf")

        for _ in range(self.n_iter):
            responsibilities, likelihood = self.expectation_step(profile)
            profile = self.maximization_step(responsibilities)
            if abs(likelihood - previous_likelihood) < self.convergence_tol:
                break
            previous_likelihood = likelihood

        responsibilities, likelihood = self.expectation_step(profile)
        starts = [int(np.argmax(gene_probs)) for gene_probs in responsibilities]
        motifs = self.extract_motifs(starts)
        return motifs, likelihood

    def choose_random_starts(self) -> List[int]:
        starts = []
        for last_index in self.last_indexes:
            starts.append(random.randint(0, last_index))
        return starts

    def extract_motifs(self, starts: List[int]) -> List[str]:
        motifs = []
        for gene, start in zip(self.genes, starts):
            motifs.append(gene[start:start + self.k])
        return motifs

    def build_profile_from_motifs(self, motifs: List[str]) -> dict:
        counts = {nucleotide: np.full(self.k, PSEUDOCOUNT, dtype=float) for nucleotide in "ACGT"}
        for motif in motifs:
            for position, nucleotide in enumerate(motif):
                counts[nucleotide][position] += 1.0
        return self.normalize_counts(counts)

    def expectation_step(self, profile: dict):
        responsibilities = []
        total_log_likelihood = 0.0

        for gene in self.genes:
            window_count = len(gene) - self.k + 1
            log_probs = np.empty(window_count, dtype=float)

            for start in range(window_count):
                log_prob = 0.0
                for position, nucleotide in enumerate(gene[start:start + self.k]):
                    log_prob += np.log(profile[nucleotide][position])
                log_probs[start] = log_prob

            max_log_prob = np.max(log_probs)
            unnormalized = np.exp(log_probs - max_log_prob)
            probs = unnormalized / np.sum(unnormalized)
            responsibilities.append(probs)
            total_log_likelihood += max_log_prob + np.log(np.sum(unnormalized))

        return responsibilities, total_log_likelihood

    def maximization_step(self, responsibilities: List[np.ndarray]) -> dict:
        counts = {nucleotide: np.full(self.k, PSEUDOCOUNT, dtype=float) for nucleotide in "ACGT"}

        for gene, gene_probs in zip(self.genes, responsibilities):
            for start, weight in enumerate(gene_probs):
                motif = gene[start:start + self.k]
                for position, nucleotide in enumerate(motif):
                    counts[nucleotide][position] += weight

        return self.normalize_counts(counts)

    @staticmethod
    def normalize_counts(counts: dict) -> dict:
        profile = {}
        total_per_column = sum(counts[nucleotide] for nucleotide in "ACGT")
        for nucleotide in "ACGT":
            profile[nucleotide] = counts[nucleotide] / total_per_column
        return profile
