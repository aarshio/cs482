import json
import random


def clean_line(line):
    allowed = set(['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
                   'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
    return "".join([c for c in line if c in allowed])


def clean(file):
    cleaned = []
    with open(file, 'r') as f:
        for line in f:
            if line[0] == '>':
                continue
            line = clean_line(line)

            while (len(line) > 40):
                if len(line) < 20:
                    break
                trunc = random.randint(20, 40)
                cleaned.append(line[:trunc])
                line = line[trunc:]

            if len(line) > 20:
                cleaned.append(line)

    return cleaned


if __name__ == '__main__':
    input_file = 'uniprot_sprot.fasta'

    kmers = {}
    allowed = set(['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
                  'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])

    print("Cleaning data...")
    data = clean(input_file)

    print("Generating kmers...")
    real_peptides = []
    fake_peptides = []
    for peptide in data:
        # real peptides
        real_peptides.append(peptide)
        for i in range(len(peptide) - 3):
            kmer = peptide[i:i + 3]
            if kmer in kmers:
                kmers[kmer][0] += 1
            else:
                kmers[kmer] = [1, 0]
        # fake peptides
        peptide_list = list(peptide)
        random.shuffle(peptide_list)
        peptide_string = "".join(peptide_list)
        fake_peptides.append(peptide_string)
        for i in range(len(peptide_string) - 3):
            kmer = peptide_string[i:i + 3]
            if kmer in kmers:
                kmers[kmer][1] += 1
            else:
                kmers[kmer] = [0, 1]

    print("Sorting kmers...")
    # sort kmers by frequency (real has higher weight), keep kmers ad dict
    kmers = sorted(
        kmers.items(), key=lambda x: x[1][0] - x[1][1], reverse=True)
    kmers = dict(kmers)

    # write data to output file
    with open("kmers.json", 'w') as outfile:
        json.dump(kmers, outfile)

    # write real peptides to real_peptides.txt
    with open("real_peptides.txt", 'w') as outfile:
        for peptide in real_peptides:
            outfile.write(peptide + "\n")

    # write fake peptides to fake_peptides.txt
    with open("fake_peptides.txt", 'w') as outfile:
        for peptide in fake_peptides:
            outfile.write(peptide + "\n")

    print("# of kmers", len(kmers.keys()))
    print("# of real and fake peptides", len(
        real_peptides), len(fake_peptides))
