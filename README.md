# fastVariant
Efficient variant caller for next generation sequencing (NGS) data

fastVariant infers variants from a bam file and outputs them in a variant call format (.vcf file).

### Input Read Strategy

### Variant Calling Strategy


### How to run?

By running ```python main.py -h``` you can see the expected input format.
```
usage: main.py [-h] --genome GENOME --sam SAM [SAM ...] [--bed BED] --rname
               RNAME --vcf VCF

Computational Biomedicine Ex 2

optional arguments:
  -h, --help           show this help message and exit
  --genome GENOME      Path to reference genome
  --sam SAM [SAM ...]  Path to SAM files. If using more than one separate
                       using a space
  --bed BED            Path to BED file
  --rname RNAME        Reference name of genome so that SAM file can be
                       filtered. Ex, use chr22 for chromosome22
  --vcf VCF            Path to VCF output file
```
An example of running command is:
```
python main.py --genome data/genome.chr22.fa --bed data/bed_chr_22.bed --sam data/alignedreads.sam --rname 'chr22' --vcf output.vcf
```
