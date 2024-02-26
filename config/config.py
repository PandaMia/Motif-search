from utils.evaluation import calc_entropy_score, calc_hamming_distance


SCORING_FUNCTIONS: dict = {
    "entropy": calc_entropy_score,
    "hamming": calc_hamming_distance,
}