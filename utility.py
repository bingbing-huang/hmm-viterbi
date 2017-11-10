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


def get_transitions():
    transitions = dict()
    transitions['Begin'] = {'State 1': 0.9999, 'State 2': 0.0001}
    transitions['State 1'] = {'State 1': 0.9999, 'State 2': 0.0001}
    transitions['State 2'] = {'State 1': 0.01, 'State 2': 0.99}
    return transitions


def main():
    exit(0)


if __name__ == '__main__':
    main()
