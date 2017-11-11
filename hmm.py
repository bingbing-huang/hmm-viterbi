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
        prob_state2_state1 = probability_matrix[1][i - 1] + np.log(transitions['State 2']['State 1']) + np.log(
            emissions['State 1'][n])
        probability_matrix[0][i] = max(prob_state1_state1, prob_state2_state1)

        prob_state1_state2 = probability_matrix[0][i - 1] + np.log(transitions['State 1']['State 2']) + np.log(
            emissions['State 2'][n])
        prob_state2_state2 = probability_matrix[1][i - 1] + np.log(transitions['State 2']['State 2']) + np.log(
            emissions['State 2'][n])
        probability_matrix[1][i] = max(prob_state1_state2, prob_state2_state2)


def perform_traceback(probability_matrix, seq, emissions, transitions, viterbi_path):
    # Traverse the probability matrix backwards one column by one column
    next_index = len(seq)
    next_state = None
    next_max_prob = None
    for prob_state1, prob_state2 in np.fliplr(probability_matrix).T:

        if next_index == len(seq):
            if prob_state1 > prob_state2:
                viterbi_path[next_index - 1] = 1
                next_state = 'State 1'
                next_max_prob = prob_state1
            else:
                viterbi_path[next_index - 1] = 2
                next_state = 'State 2'
                next_max_prob = prob_state2
            next_index -= 1
            continue

        if prob_state1 + np.log(transitions['State 1'][next_state]) + np.log(
                emissions[next_state][seq[next_index]]) == next_max_prob:
            viterbi_path[next_index - 1] = 1
            next_state = 'State 1'
            next_max_prob = prob_state1
        else:
            viterbi_path[next_index - 1] = 2
            next_state = 'State 2'
            next_max_prob = prob_state2
        next_index -= 1


def find_hits(viterbi_path):
    # 0 based index
    hit_indices = list()
    hit_found = False

    hit_start_index = None
    for index, state in enumerate(viterbi_path):

        if state == 2 and not hit_found:
            hit_start_index = index
            hit_found = True

        if hit_found and state == 1:
            hit_indices.append((hit_start_index, index - 1))
            hit_found = False

    return hit_indices


def main():
    seq = load_sequence()
    emissions = get_emissions()
    transitions = get_transitions()
    """
        state 1 (row 0): ...
        state 2 (row 1: ...
    """
    probability_matrix = np.empty((2, len(seq)), dtype='float64')
    calculate_probabilities(seq, emissions, transitions, probability_matrix)
    viterbi_path = np.empty(len(seq), dtype='int')
    perform_traceback(probability_matrix, seq, emissions, transitions, viterbi_path)
    hit_indices = find_hits(viterbi_path)


if __name__ == '__main__':
    main()
