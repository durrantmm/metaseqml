import sys
import gzip
import itertools
import re
from os.path import dirname, basename, join, splitext
from snakemake import shell

def main(fastq, taxonomy_path, outfile, k):

    shell('mkdir -p {dir}'.format(dir=dirname(outfile)))

    bases = ['A', 'T', 'G', 'C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]



    taxon = basename(fastq).split('.')[0]
    hierarchy = get_hierarchy_names(taxonomy_path, taxon)

    with gzip.open(fastq) as filein:

        with open(splitext(outfile)[0], 'w') as fileout:

            fileout.write("\t".join(["Seq"] + ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'] + kmers) + '\n')

            count = 0
            for index, line in enumerate(filein):
                if count == 10000:
                    break
                line = line.decode("UTF-8")

                if index % 4 == 1:

                    seq = line.strip()

                    kmer_comp_dict = get_kmer_comp(seq, k, kmers)
                    kmer_comp = str_kmer_comp(seq, kmer_comp_dict, kmers)
                    count += 1
                    fileout.write("\t".join([seq] + hierarchy + kmer_comp)+'\n')

    shell('gzip {out}'.format(out=splitext(outfile)[0]))



def get_hierarchy_names(taxonomy_path, taxon):
    nodes_path = join(taxonomy_path, 'nodes.dmp')
    merged_path = join(taxonomy_path, 'merged.dmp')
    names_path = join(taxonomy_path, 'names.dmp')

    nodes, ranks = get_taxon_nodes_ranks(nodes_path, merged_path)
    names = taxa_to_names(names_path)

    hierarchy = get_taxon_hierarchy_list(taxon, nodes)
    hierarchy_names = taxon_list_to_names(hierarchy, ranks, names)
    return hierarchy_names


def taxa_to_names(names_path):

    taxa2names = dict()
    with open(names_path) as infile:
        for line in infile:
            taxon_id, name, junk, name_type = line.strip('\t|\n').split('\t|\t')
            if name_type == 'scientific name':
                taxa2names[taxon_id] = name
    return taxa2names


def taxon_list_to_names(taxon_list, ranks, names):
    keep_ranks = [['superkingdom', 'kingdom'],
                  ['phylum'],
                  ['class'],
                  ['order'],
                  ['family'],
                  ['genus'],
                  ['species']]
    outlist = []
    existing_ranks = []
    for rank in keep_ranks:
        i=0
        while i < len(taxon_list) and ranks[taxon_list[i]] not in rank:
            i += 1
        if i == len(taxon_list):
            outlist.append(outlist[-1])
        else:
            outlist.append(names[taxon_list[i]])

        existing_ranks.append(rank)

    return outlist


def get_taxon_hierarchy_list(taxon_id, taxon_nodes_dict):
    hierarchy = [taxon_id]

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.append(taxon_id)

    return list(reversed(hierarchy))


def get_taxon_nodes_ranks(nodes_path, merged_path):

    taxon_nodes_dict = {}
    taxon_ranks_dict = {}
    with open(nodes_path) as nodes_in:
        for line in nodes_in:
            line = line.strip().split("|")
            id = line[0].strip()
            rank = line[2].strip()
            parent_id = line[1].strip()
            taxon_nodes_dict[id] = parent_id
            taxon_ranks_dict[id] = rank

    with open(merged_path) as merged_in:
        for line in merged_in:
            line = line.strip().split("|")
            orig_id = line[0].strip()
            merged_id = line[1].strip()
            taxon_nodes_dict[orig_id] = merged_id
            taxon_ranks_dict[orig_id] = taxon_ranks_dict[merged_id]

    return taxon_nodes_dict, taxon_ranks_dict



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
    return list(map(str, [kmer_dict[kmer] for kmer in kmers]))

if __name__ == '__main__':

    fastq_path = sys.argv[1]
    taxonomy_path = sys.argv[2]
    outfile = sys.argv[3]
    k = int(sys.argv[4])
    main(fastq_path, taxonomy_path, outfile, k)