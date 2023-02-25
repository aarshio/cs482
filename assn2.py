import argparse
import json
import math
import numpy as np
from sklearn.metrics import roc_auc_score


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


def score1(peptide, freq_table):
    sc1 = 0
    for i in range(len(peptide) - 3):
        kmer = peptide[i:i + 3]
        sc1 += math.log(freq_table[kmer][0] / freq_table[kmer][1], 2)
    return sc1


# evaluate AUROC score
def evaluate(peptides, freq_table):
    # get real_peptides
    real_peptides = []
    with open('real_peptides.txt', 'r') as f:
        for line in f:
            real_peptides.append(line.strip())

    # get fake_peptides
    fake_peptides = []
    with open('fake_peptides.txt', 'r') as f:
        for line in f:
            fake_peptides.append(line.strip())

    # assign labels to real_peptides
    peptides = []
    for peptide in real_peptides:
        sc1 = score1(peptide, freq_table)
        p1 = 1 / (1 + 2**(-sc1))
        print(sc1, '\t', p1)
        peptides.append((peptide, 1, p1))
    # assign labels to fake_peptides
    for peptide in fake_peptides:
        sc1 = score1(peptide, freq_table)
        p1 = 1 / (1 + 2**(-sc1))
        print(sc1, '\t', p1)
        peptides.append((peptide, 0, p1))

    # calculate AUROC score
    y_true = [p[1] for p in peptides]
    y_score = [p[2] for p in peptides]
    return roc_auc_score(y_true, y_score)


if __name__ == '__main__':
    args = read_arguments()
    peptides = read_input_file(args.input)
    # load frequency table from cache.json
    freq_table = {}
    with open('kmers.json', 'r') as f:
        freq_table = json.load(f)

    ans = []

    observed = {}
    for peptide in peptides:
        ans.append((peptide, score1(peptide, freq_table)))

    for a in ans:
        print(a[1], '\t', a[0])

    # print the AUROC score
    print(evaluate(peptides, freq_table))
