import reader as io
import utils
from variant_caller import naive_caller
import argparse

# TODO:
# Add SVLEN, SVTYPE information to output -- DONE
# Make deletions work (insertions work fine) -- DONE
# Accept multiple input files --DONE by Tommaso
# Make processing of reads faster when a BED file is given
# Some rsXXXXXXXX information is found incorrectly from the BED file
# # Possibly search using maximum overlap instead of simple dict? (slower)

# python main.py --genome genome.chr22.fa --sam alignedreads.sam --rname chr22 --vcf output.vcf
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Computational Biomedicine Ex 2')
    parser.add_argument(
        '--genome', help='Path to reference genome', required=True)
    parser.add_argument('--sam', nargs="+", help='Path to SAM files. If using more than one separate using a space', required=True)
    parser.add_argument('--bed', help='Path to BED file', required=False)
    parser.add_argument(
        '--rname',
        help='Reference name of genome so that SAM file can be filtered. Ex, use chr22 for chromosome22',
        required=True)
    parser.add_argument('--vcf', help='Path to VCF output file', required=True)

    args = parser.parse_args()

    reader = io.IO(args.sam, args.rname, args.genome, args.vcf, args.bed)
    reads = reader.return_reads()
    bed = reader.read_bed()
    description, short_name, genome = reader.ref_seq()

    candidates = utils.countAllReads(genome, reads)

    # utils.pickle_result(bed, "b_file")
    # utils.pickle_result(genome, "g_file")
    # utils.pickle_result(candidates, "candidates_file")
    # candidates = utils.open_jar("candidates_file")
    # genome = utils.open_jar("g_file")
    # bed = utils.open_jar("b_file")

    # variant caller
    # print(candidates)
    vc_list = naive_caller(candidates, genome, bed)
    # output
    reader.produce_to_file(vc_list, short_name, 'SVLEN=1')
