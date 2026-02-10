from typing import List
import random
import numpy as np
from tqdm import tqdm
from motif_search.base_search import BaseMotifSearch
from utils.motif_utils import find_consensus
from utils.decorators import error_handler


class GeneticMotifSearch(BaseMotifSearch):
    """Genetic algorithm to find best Motif in DNA sequences"""

    def __init__(
        self, 
        genes: List[str], 
        k: int, 
        metric: str, 
        n_iter: int, 
        n_bots: int, 
        n_surv: int, 
        mutation_coef: float
    ):
        """Initialize genetic algorithm for motif search.

        Args:
            genes (list): list of genes: genes[str]
            k (int): length of the motif
            metric (str): metric to evaluate found motifs
            n_iter (int): number of steps of genetic.
            n_bots (int): population size on each epoch.
            n_surv (int): number of survivors after epoch.
            mutation_coef (float): coefficient of mutation for each position in new bot.
        """
        super().__init__(genes, k, metric)        # initialize BaseMotifSearch
        self.n_iter = n_iter                      # number of epochs
        self.n_bots = n_bots                      # population size
        self.n_surv = n_surv                      # number of survivors
        self.n_new_bots = n_bots - n_surv         # number of new bots
        self.mutation_coef = mutation_coef        # mutation coefficient
        self.last_index = self.gene_len - self.k  # last possible index for motif start

    error_handler
    def run_search(self) -> dict:
        """Run genetic algorithm.

        Returns:
            dict: response with found motifs for each gene, their scores and consensus motif
        """
        # run genetic
        population = self.create_first_population()
        scores = self.evaluate_population(population)
        best_bot = None
        best_score = float("inf")
        for _ in tqdm(range(self.n_iter)):
            min_idx, min_score = min(enumerate(scores), key=lambda pair: pair[1])
            if min_score < best_score:
                best_score = min_score
                best_bot = population[min_idx]
            population = self.create_next_generation(population, scores)
            scores = self.evaluate_population(population)
        if best_bot is None:
            min_idx, min_score = min(enumerate(scores), key=lambda pair: pair[1])
            best_score = min_score
            best_bot = population[min_idx]

        # choose best motifs
        best_motifs = [self.genes[i][start_idx:start_idx+self.k] for i, start_idx in enumerate(best_bot)]

        # get scores for the best motifs and consensus motif
        scores = self.evaluate_best_motifs(best_motifs)
        consensus = find_consensus(best_motifs)
        
        response = {
            "best_motifs": best_motifs,
            "scores": scores,
            "consensus": consensus
        }
        
        return response
    
    def create_random_bot(self):
        bot = [None] * self.n_genes
        for i in range(self.n_genes):
            bot[i] = random.randint(0, self.last_index)
        return bot
    
    def create_first_population(self):
        population = []
        for _ in range(self.n_bots):
            bot = self.create_random_bot()
            population.append(bot)
        return population
    
    def choose_survivors(self, population, scores):
        sorted_population = [bot for _, bot in sorted(zip(scores, population), key=lambda pair: pair[0])]
        survivors = sorted_population[:self.n_surv]
        return survivors

    def create_next_generation(self, old_population, scores):
        survivors = self.choose_survivors(old_population, scores)
        new_population = []
        for _ in range(self.n_new_bots):
            parents_ids = np.random.choice(self.n_surv, 2, replace=False).tolist()
            parent_1 = survivors[parents_ids[0]]
            parent_2 = survivors[parents_ids[1]]
            new_bot = [None] * self.n_genes
            for i in range(self.n_genes):
                if random.random() < self.mutation_coef:
                    new_bot[i] = random.randint(0, self.last_index)
                elif random.random() < 0.5:
                    new_bot[i] = parent_1[i]
                else:
                    new_bot[i] = parent_2[i]
            new_population.append(new_bot)
        return new_population

    def evaluate_population(self, population):
        scores = []
        for bot in population:
            motifs = [self.genes[gene_idx][start_idx:start_idx+self.k] for gene_idx, start_idx in enumerate(bot)]
            score = self.scoring_function(motifs)
            scores.append(score)
        return scores
