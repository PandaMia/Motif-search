from typing import List, Optional, Literal, Dict, Any
from pydantic import BaseModel


class FindMotif(BaseModel):
    genes: List[str]
    k: int
    method: Optional[Literal["greedy", "random", "gibbs", "genetic"]] = "gibbs"
    metric: Optional[Literal["hamming", "entropy"]] = "hamming"
    n_iter: Optional[int] = 1000
    n_bots: Optional[int] = 100
    n_surv: Optional[int] = 20
    mutation_coef: Optional[float] = 0.15


class FindMotifResponse(BaseModel):
    best_motifs: Optional[List[str]] = None
    scores: Optional[Dict[str, Any]] = None
    consensus: Optional[str] = None
    error: Optional[str] = None
