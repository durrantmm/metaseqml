import os, sys
from glob import glob
from os.path import join, dirname, basename
from shutil import copyfile

def main(midas_rep_genomes_path, taxonomy_path):

    genomes = glob(join(midas_rep_genomes_path, '**','genome.fna.gz'))
    names2taxa = names_to_taxon_ids(taxonomy_path)

    taxa2paths = taxa_to_paths(genomes, names2taxa)
    out_dir = join(dirname(midas_rep_genomes_path), '0.genomes')
    os.makedirs(out_dir)

    for taxon in taxa2paths:
        path = taxa2paths[taxon]
        copyfile(path, join(out_dir, '{taxon}.fasta.gz'.format(taxon=taxon)))



def taxa_to_paths(genomes, names2taxa):
    taxa2paths = dict()
    for path in genomes:
        name = ' '.join(path.split('/')[-2].split('_')[:-1])
        try:
            taxa2paths[names2taxa[name]] = path
        except:
            continue
    return taxa2paths


def names_to_taxon_ids(taxonomy_path):

    names_path = join(taxonomy_path, 'names.dmp')

    names2taxa = dict()
    with open(names_path) as infile:
        for line in infile:
            taxon_id, name, junk, name_type = line.strip('\t|\n').split('\t|\t')
            if name_type == 'scientific name':
                names2taxa[name] = taxon_id
    return names2taxa



if __name__ == '__main__':
    midas_rep_genomes_path = sys.argv[1]
    taxonomy_path = sys.argv[2]

    main(midas_rep_genomes_path, taxonomy_path)
