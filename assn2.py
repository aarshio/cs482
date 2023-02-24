import argparse
import json
import math


def read_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input file', required=True)
    parser.add_argument('--output', help='output file', required=True)
    args = parser.parse_args()
    return args


def read_input_file(input_file):
    with open(input_file, 'r') as f:
        seqs = []
        for line in f:
            seqs.append(line.strip())
        return seqs


if __name__ == '__main__':
    args = read_arguments()
    peptides = read_input_file(args.input)
    # load frequency table from cache.json
    real_table = {}
    with open('real.json', 'r') as f:
        real_table = json.load(f)

    fake_table = {}
    with open('fake.json', 'r') as f:
        fake_table = json.load(f)

    ans = []

    observed = {}
    for peptide in peptides:
        score1 = 0
        for i in range(len(peptide) - 3):
            kmer = peptide[i:i + 3]

            score1 += math.log(real_table[kmer] / fake_table[kmer], 2)

        ans.append((peptide, score1))

    for a in ans:
        print(a[0], a[1])
