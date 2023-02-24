import json
import random

if __name__ == '__main__':
    input_file = 'uniprot_sprot.fasta'

    real_peptides = []
    kmers = {}
    allowed = set(['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])

    with open(input_file, 'r') as f:
        peptide = ''
        for line in f:
            if line[0] == '>':
                if peptide == '':
                    continue

                # # randomly truncate peptides to length between 20 and 40
                # length = random.randint(20, 40)
                # peptide = peptide[:length]

                for i in range(len(peptide) - 3):
                    kmer = peptide[i:i + 3]
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
                real_peptides.append(peptide)
                peptide = ''
            else:
                # replace any letters with "" if not in allowed
                line = "".join([c for c in line if c in allowed])
                peptide += line

    # sort kmers by value
    kmers = {k: v for k, v in sorted(kmers.items(), key=lambda item: item[1], reverse=True)}

    # # normalize kmers
    # total = sum(kmers.values())
    # for kmer in kmers:
    #     kmers[kmer] /= total

    # write data to output file
    with open("real.json", 'w') as outfile:
        json.dump(kmers, outfile)

    # write real peptides to real_peptides.txt
    with open("real_peptides.txt", 'w') as outfile:
        for peptide in real_peptides:
            outfile.write(peptide + "\n")

    print(len(kmers.keys()))
    print(len(real_peptides))
