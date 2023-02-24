import json
import random

# sample half of data from real data (uniprot_sprot.fasta)
# and create other half peptides uniformly randomly distributed between 20 and 40 (inclusive)
# put kmer (3) frequency count into fake.json


def read_input_file(input_file):
    with open(input_file, 'r') as f:
        seqs = []
        for line in f:
            if line[0] == '>':
                continue
            seqs.append(line.strip())
        return seqs


if __name__ == '__main__':
    input_file = 'real_peptides.txt'

    # get all peptides from input file
    real_peptides = read_input_file(input_file)

    allowed = set(['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])

    sample = []
    for peptide in real_peptides:
        # append shuffled peptides to sample
        peptide = "".join([c for c in peptide if c in allowed])
        string_list = list(peptide)
        random.shuffle(string_list)
        shuffled_string = "".join(string_list)
        sample.append(shuffled_string)

    kmers = {}

    for peptide in sample:
        # replace any letters with "" if not in allowed
        peptide = "".join([c for c in peptide if c in allowed])
        for i in range(len(peptide) - 3):
            kmer = peptide[i:i + 3]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1

    # sort kmers by value
    kmers = {k: v for k, v in sorted(kmers.items(), key=lambda item: item[1], reverse=True)}

    # # normalize kmers
    # total = sum(kmers.values())
    # for kmer in kmers:
    #     kmers[kmer] /= total

    # write data to output file
    with open("fake.json", 'w') as outfile:
        json.dump(kmers, outfile)

    # write sample to fake_peptides.txt
    with open("fake_peptides.txt", 'w') as outfile:
        for peptide in sample:
            outfile.write(peptide + "\n")

    print(len(kmers.keys()))
    print(len(sample))
