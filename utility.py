from Bio import SeqIO
import re


def get_sequence():
    fasta_sequences = SeqIO.parse(open('GCF_000091665.1_ASM9166v1_genomic.fna', 'r'), 'fasta')
    for sequence in fasta_sequences:
        sequence = str(sequence.seq)
        sequence.replace(" ", "")
        sequence.upper()
        return re.sub(r'([^ACGT])', "T", sequence)


def get_cds():
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


def main():
    exit(0)


if __name__ == '__main__':
    main()