# Motif Search Methods

Short descriptions of the motif-search methods available in the service.

## Greedy

Builds motifs step by step from the most promising local choice.
It is fast and deterministic, but it can get stuck in a locally good solution instead of the global best motif.

## Random

Starts from random motif positions and repeatedly improves them with the current profile matrix.
It is simple and often effective, but the final answer depends on the starting points, so more iterations usually improve stability.

## Gibbs

Uses Gibbs sampling to update one sequence at a time while keeping the others fixed.
This stochastic method explores the search space better than purely greedy updates and can escape some local optima.

## Genetic

Represents motif positions as a population of candidate solutions and evolves them through selection, crossover, and mutation.
It is useful for broad search over difficult datasets, but it usually needs more tuning and compute than simpler methods.

## EM

Uses expectation-maximization with a probabilistic motif profile and multiple random restarts.
It is a good fit when you want a profile-based method that is more principled than heuristic local search, but it can still converge to a local optimum without enough restarts.
