from Bio import SeqIO
import re


def load_sequence():
    fasta_sequences = SeqIO.parse(open('GCF_000091665.1_ASM9166v1_genomic.fna', 'r'), 'fasta')
    for sequence in fasta_sequences:
        sequence = str(sequence.seq)
        sequence.replace(" ", "")
        sequence.upper()
        return re.sub(r'([^ACGT])', "T", sequence)


def load_cds():
    gen_bank = SeqIO.parse(open('GCF_000091665.1_ASM9166v1_genomic.gbff', 'r'), 'genbank')
    for records in gen_bank:
        ignored_chars = ['>', '<']
        end_indices = set()
        for feature in records.features:
            if feature.type == 'CDS' and feature.strand == 1:
                start = str(feature.location.start)
                end = str(feature.location.end)
                if any(x in start for x in ignored_chars) or any(x in end for x in ignored_chars):
                    continue
                # convert 0 based index to 1 based index for start index
                end_indices.add(int(end))
        return end_indices


def get_emissions():
    emissions = dict()
    emissions['State 1'] = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    emissions['State 2'] = {'A': 0.20, 'C': 0.30, 'G': 0.30, 'T': 0.20}
    return emissions


def print_emission(emissions):
    print("The HMM emission parameters:")
    line = '{:13} {:5} {:5} {:5} {:5}'.format('Emissions', 'A', 'C', 'G', 'T')
    print(line)
    line = '{:12} {:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format('State 1', emissions['State 1']['A'],
                                                           emissions['State 1']['C'],
                                                           emissions['State 1']['G'], emissions['State 1']['T'])
    print(line)
    line = '{:12} {:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format('State 2', emissions['State 2']['A'],
                                                           emissions['State 2']['C'],
                                                           emissions['State 2']['G'], emissions['State 2']['T'])
    print(line)


def get_transitions():
    transitions = dict()
    transitions['Begin'] = {'State 1': 0.9999, 'State 2': 0.0001}
    transitions['State 1'] = {'State 1': 0.9999, 'State 2': 0.0001}
    transitions['State 2'] = {'State 1': 0.01, 'State 2': 0.99}
    return transitions


def print_transitions(transitions):
    print("The HMM transitions parameters:")
    line = '{:15}  {:11}  {:11}'.format('Transitions', 'State 1', 'State 2')
    print(line)
    line = '{:11}  {:13.7f}  {:11.7f}'.format('Begin', transitions['Begin']['State 1'], transitions['Begin']['State 2'])
    print(line)
    line = '{:11}  {:13.7f}  {:11.7f}'.format('State 1', transitions['State 1']['State 1'],
                                           transitions['State 1']['State 2'])
    print(line)
    line = '{:11}  {:13.7f}  {:11.7f}'.format('State 2', transitions['State 2']['State 1'],
                                           transitions['State 2']['State 2'])
    print(line)


def print_log_probability_viterbi(probability_matrix):
    _, columns = probability_matrix.shape
    print("The log probability of the overall Viterbi path: %f" % max(probability_matrix[0][columns - 1],
                                                                     probability_matrix[1][columns - 1]))


def print_hits(hit_indices):
    print("The total number of hits found: %d" % len(hit_indices))


def print_hits_details(number_of_hits_to_print, length_hits, hit_indices):
    sorted_hit_indices = sorted(hit_indices, key=lambda x: x[0])
    for hit_start, hit_end in sorted_hit_indices:
        if number_of_hits_to_print == 0:
            break
        length = hit_end - hit_start + 1
        if length >= length_hits:
            print("Start: %d, End: %d, Length: %d" % (hit_start + 1, hit_end + 1, length))
            number_of_hits_to_print -= 1


def print_language_bake_off(process_times, hit_indices):
    print("Here is the language bake-off:")
    print("\nPython 3.6.2")
    print("\nThe total CPU time taken by the algorithm for the first 9 Viterbi iterations: %fs" % sum(process_times[:len(process_times) - 1]))
    print("\nThe basic info about the computer: 2.4 GHz Intel Core i7, 8 GB RAM")
    print("\nThe first 10 hits of length >= 50 from the final iteration:")
    print_hits_details(10, 50, hit_indices)


def main():
    exit(0)


if __name__ == '__main__':
    main()
