import sys
from os.path import splitext
from snakemake import shell

FASTQ_DIR = sys.argv[1]
input = sys.argv[2]
output_f1 = sys.argv[3]
output_f2 = sys.argv[4]

output_f1_noext = splitext(output_f1)[0]
output_f2_noext = splitext(output_f2)[0]
shell("mkdir -p {fastq_dir}".format(fastq_dir=FASTQ_DIR))
shell("wgsim -1 100 -2 100 {input} {output_f1} {output_f2}".format(input=input, output_f1=output_f1_noext, output_f2=output_f2_noext))
shell("gzip {output_f1} {output_f2}".format(output_f1=output_f1_noext, output_f2=output_f2_noext))