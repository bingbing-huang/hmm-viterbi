from hmm_viterbi import *
from utility import *
import time


def main():
    seq = load_sequence()
    emissions = get_emissions()
    transitions = get_transitions()
    iteration_times = 10
    process_times = []
    hit_indices = None

    for i in range(iteration_times):
        """
            state 1 (row 0): ...
            state 2 (row 1: ...
        """
        probability_matrix = np.empty((2, len(seq)), dtype='float64')
        viterbi_path = np.empty(len(seq), dtype='int')

        start_time = time.process_time()
        n_count = calculate_probabilities(seq, emissions, transitions, probability_matrix)
        perform_traceback(probability_matrix, seq, emissions, transitions, viterbi_path)
        hit_indices = find_hits(viterbi_path)
        process_times.append(time.process_time() - start_time)

        print("Iteration %d:" % (i + 1))
        print()

        print_emission(emissions)
        print()

        print_transitions(transitions)
        print()

        print_log_probability_viterbi(probability_matrix)
        print()

        print_hits(hit_indices)
        print()

        print("Here is detail information for the hits: ")
        if i == 9:
            print_hits_details(len(hit_indices), 0, hit_indices)
        else:
            print_hits_details(5, 0, hit_indices)

        print()
        print("***********************************************************")
        print()
        update_emissions(emissions, seq, n_count, hit_indices)
        update_transitions(transitions, n_count, hit_indices)

    print_language_bake_off(process_times, hit_indices)


if __name__ == '__main__':
    main()
