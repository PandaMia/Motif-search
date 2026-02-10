from config import base_models
from utils.decorators import error_handler
from motif_search.greedy_search import GreedyMotifSearch
from motif_search.random_search import RandomMotifSearch
from motif_search.genetic_search import GeneticMotifSearch
from motif_search.gibbs_search import GibbsMotifSearch


class Manager():
    def __init__(self, request_data: base_models.FindMotif):
        self.genes = request_data.genes
        self.k = request_data.k
        self.metric = request_data.metric
        self.params = request_data

    @error_handler
    def find_motif(self):
        """Find best motif in DNA sequences.

        Returns:
            Tuple[List[str], dict]: found best motifs and scores
        """
        if self.params.method == "greedy":
            self.prepare_greedy()
        elif self.params.method == "random":
            self.prepare_random()
        elif self.params.method == "gibbs":
            self.prepare_gibbs()
        elif self.params.method == "genetic":
            self.prepare_genetic()
        response = self.search.run_search()
        return response
    
    def prepare_greedy(self):
        self.search = GreedyMotifSearch(
            self.genes, 
            self.k, 
            self.metric
        )

    def prepare_random(self):
        n_iter = self.params.n_iter
        self.search = RandomMotifSearch(
            self.genes, 
            self.k,  
            self.metric,
            n_iter
        )

    def prepare_gibbs(self):
        n_iter = self.params.n_iter
        self.search = GibbsMotifSearch(
            self.genes, 
            self.k, 
            self.metric, 
            n_iter
        )
    
    def prepare_genetic(self):
        n_iter = self.params.n_iter 
        n_bots = self.params.n_bots 
        n_surv = self.params.n_surv 
        mutation_coef = self.params.mutation_coef
        self.search = GeneticMotifSearch(
            self.genes, 
            self.k, 
            self.metric, 
            n_iter, 
            n_bots, 
            n_surv, 
            mutation_coef
        )
