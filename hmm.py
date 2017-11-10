from utility import *
import numpy as np


def calculate_probabilities(seq, emissions, transitions, probability_matrix):
    """
        state 1 (row 0): ...
        state 2 (row 1: ...
    """
    for i, n in enumerate(seq):
        # Initialize probabilities for the first nucleotide
        if i == 0:
            probability_matrix[0][i] = np.log(transitions['Begin']['State 1']) + np.log(emissions['State 1'][n])
            probability_matrix[1][i] = np.log(transitions['Begin']['State 2']) + np.log(emissions['State 2'][n])
            continue

        prob_state1_state1 = probability_matrix[0][i - 1] + np.log(transitions['State 1']['State 1']) + np.log(
            emissions['State 1'][n])
        prob_state2_state1 = probability_matrix[1][i - 1] + np.log(transitions['State 2']['State 2']) + np.log(
            emissions['State 1'][n])
        probability_matrix[0][i] = max(prob_state1_state1, prob_state2_state1)

        prob_state1_state2 = probability_matrix[0][i - 1] + np.log(transitions['State 1']['State 2']) + np.log(
            emissions['State 2'][n])
        prob_state2_state2 = probability_matrix[1][i - 1] + np.log(transitions['State 2']['State 2']) + np.log(
            emissions['State 2'][n])
        probability_matrix[1][i] = max(prob_state1_state2, prob_state2_state2)


def main():
    seq = load_sequence()
    emissions = get_emissions()
    transitions = get_transitions()
    """
        state 1 (row 0): ...
        state 2 (row 1: ...
    """
    probability_matrix = np.zeros((2, len(seq)), dtype='float64')
    calculate_probabilities(seq, emissions, transitions, probability_matrix)
    print(probability_matrix)


if __name__ == '__main__':
    main()
