import sys
import gzip
import itertools
import re

def main(fastq):
    bases = ['A', 'T', 'G', 'C']
    k = 3
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]

    print("\t".join(["Seq"] + kmers))

    with gzip.open(fastq) as filein:

        for index, line in enumerate(filein):
            if index % 4 == 1:

                seq = line.strip()
                kmer_comp_dict = get_kmer_comp(seq, k, kmers)
                kmer_comp = str_kmer_comp(seq, kmer_comp_dict, kmers)
                print("\t".join([seq] + kmer_comp))



def get_kmer_comp(seq, k, kmers):
    seq = trim_seq_N(seq)

    kmer_dict = {kmer: 0 for kmer in kmers}

    for i in range(0, len(seq)-k):
        kmer_dict[seq[i:i+k]] += 1
    n_kmers = float(len(seq)-k)
    kmer_dict = {kmer: kmer_dict[kmer]/n_kmers for kmer in kmer_dict}

    return kmer_dict

def trim_seq_N(seq):

    if seq.find('N') > -1:
        pos = [n.start() for n in re.finditer('N', seq)]

        if 0 not in pos:
            pos.insert(0, 0)

        if len(seq)not in pos:
            pos.append(len(seq))

        largest_seq_length = 0
        largest_seq = seq
        for p in range(len(pos)-1):

            length = pos[p+1] - pos[p]
            if length > largest_seq_length:
                largest_seq_length = length
                largest_seq = seq[pos[p]+1: pos[p+1]]
        seq = largest_seq

    return seq

def str_kmer_comp(seq, kmer_dict, kmers):
    return map(str, [kmer_dict[kmer] for kmer in kmers])

if __name__ == '__main__':

    fastq_path = sys.argv[1]
    main(fastq_path)