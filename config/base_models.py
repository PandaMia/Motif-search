from typing import List, Optional
from pydantic import BaseModel


class FindMofif(BaseModel):
    genes: List[str]
    k: int
    method: Optional[str] = "gibbs"
    metric: Optional[str] = "hamming"
    n_iter: Optional[int] = 1000
    n_bots: Optional[int] = 100
    n_surv: Optional[int] = 20
    mutation_coef: Optional[float] = 0.15
