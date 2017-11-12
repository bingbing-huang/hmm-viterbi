from utility import *
import numpy as np


def calculate_probabilities(seq, emissions, transitions, probability_matrix):
    """
        state 1 (row 0): ...
        state 2 (row 1: ...
    """
    # count of nucleonide
    n_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for i, n in enumerate(seq):
        n_count[n] += 1
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
    return n_count


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
        # The end index is the last index of the continuous state 2 nucleotides
        if hit_found and state == 1:
            hit_indices.append((hit_start_index, index - 1))
            hit_found = False

    return hit_indices


def update_emissions(emissions, seq, n_count, hit_indices):
    # update for state 2
    state2_n_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for hit_start, hit_end in hit_indices:
        for n in seq[hit_start:hit_end + 1]:
            state2_n_count[n] += 1
    state2_total_count = sum(state2_n_count.values())
    for n, count in state2_n_count.items():
        emissions['State 2'][n] = count / state2_total_count

    # update for state 1
    state1_total_count = sum(n_count.values()) - state2_total_count
    for n, count in n_count.items():
        emissions['State 1'][n] = (count - state2_n_count[n]) / state1_total_count


def update_transitions(transitions, n_count, hit_indices):
    total_number_transitions = sum(n_count.values()) - 1

    count_state1_state2 = 0
    count_state2_state1 = 0
    count_state2_state2 = 0
    for hit_start, hit_end in hit_indices:

        if hit_start != 0:
            count_state1_state2 += 1

        if hit_end != total_number_transitions:
            count_state2_state1 += 1

        count_state2_state2 += hit_end - hit_start
    count_state1_state1 = total_number_transitions - (count_state1_state2 + count_state2_state1 + count_state2_state2)

    count_state1 = count_state1_state1 + count_state1_state2
    transitions['State 1']['State 1'] = count_state1_state1 / count_state1
    transitions['State 1']['State 2'] = count_state1_state2 / count_state1
    count_state2 = count_state2_state1 + count_state2_state2
    transitions['State 2']['State 1'] = count_state2_state1 / count_state2
    transitions['State 2']['State 2'] = count_state2_state2 / count_state2


def main():
    seq = load_sequence()
    emissions = get_emissions()
    transitions = get_transitions()
    iteration_times = 10
    """
        state 1 (row 0): ...
        state 2 (row 1: ...
    """
    for i in range(iteration_times):
        probability_matrix = np.empty((2, len(seq)), dtype='float64')
        n_count = calculate_probabilities(seq, emissions, transitions, probability_matrix)
        viterbi_path = np.empty(len(seq), dtype='int')
        perform_traceback(probability_matrix, seq, emissions, transitions, viterbi_path)
        hit_indices = find_hits(viterbi_path)

        print("Iteration %d:" % (i + 1))
        print()

        update_emissions(emissions, seq, n_count, hit_indices)
        print_emission(emissions)
        print()

        update_transitions(transitions, n_count, hit_indices)
        print_transitions(transitions)
        print()

        print_log_probability_viterbi(probability_matrix)
        print()

        print_hits(hit_indices)
        print()

        if i == 9:
            print_hits_details(len(hit_indices), 0, hit_indices)
        else:
            print_hits_details(5, 0, hit_indices)

        print()
        print("***********************************************************")
        print()


if __name__ == '__main__':
    main()
