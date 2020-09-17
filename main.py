"""Console script for converting a VCF into a gVCF."""
import argparse
import logging
import sys
import utils.bed_reader as br
import utils.vcf_mangler as v
from pyfaidx import Fasta


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcf", dest='vcf', required=True, help="VCF File to convert")
    parser.add_argument("-b", "--bed", dest='bed', help="Target BED file", required=True)
    parser.add_argument("-f", "--fasta", dest='fasta', required=True,
                        help="Where is the reference genome (FASTA) file?")
    parser.add_argument("-o", "--outname", dest='outname', default='output.g.vcf',
                        help="What should the output file name be? [output.g.vcf]")
    parser.add_argument("-M", "--min_depth", dest='min_depth', default='20', help="Dummy value for min_depth [20]")
    parser.add_argument("-P", "--pl", dest='pl', default='99', help="Dummy value for PL (Non zero genotypes) [99]")
    parser.add_argument("-V", "--verbose", dest="logLevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default="INFO", help="Set the logging level [INFO]")
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stderr, level=args.logLevel, format='%(name)s (%(levelname)s): %(message)s')
    logger = logging.getLogger(__name__)
    logger.setLevel(args.logLevel)
    return args


if __name__ == "__main__":
    args = parse_args()
    logging.debug(f'{args}')
    # Read the the BED file, and convert its input into a list of lists
    bed_regions = br.parse_bed(args.bed)
    # Create the FASTA Object for variant extraction
    fasta = Fasta(args.fasta)
    history = dict()
    history['current'] = None
    outfile = open(args.outname, 'w')
    with open(args.vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith("#CHROM"):
                    outfile.write('##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within '
                                  'the GVCF block">' + '\n')
                    outfile.write('##INFO=<ID=Inferred,Number=1,Type=Bool,Description="Whether the gVCF value was '
                                  'inferred from a VCF and BED file rather than a true gVCF">' + '\n')
                    outfile.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">'
                                  + '\n')
                    assert len(line.strip().split('\t')) < 11, "This script only works on single sample VCFs!"
                    outfile.write(line)
                else:
                    # Regular header lines
                    outfile.write(line)
            else:
                # Parsing individual records
                line_elements = line.strip().split('\t')

                # Check for multiallele and just write out
                if line_elements[4].__contains__(','):
                    outfile.write(line + '\n')

                # Check to see if this is the first time
                if history['current'] is None:
                    history['current'] = line_elements[0]
                    a, b, c, d = br.get_chrom_regions(history['current'], bed_regions)
                    history['current_max'] = int(a)
                    history['current_min'] = int(b)
                    history['current_regions'] = c

                # If its not the first variant, but is the same chromosome,
                if history['current'] == line_elements[0]:
                    # Maintain an account of min and max regions for that chrom
                    # If variant pos is < max region,
                    #   Create a region from the min to the variant with pre-specified values
                    #   Add the Variant and the NON-REF (with appropriate FORMAT FIELDS)
                    #   Reset chrom minimum
                    current_pos = int(line_elements[1])
                    assert current_pos > history['current_min'], "current_pos < history['current_min']. " \
                                                                 "Is your VCF sorted?"
                    if current_pos < history['current_max']:
                        logging.warning(
                            f"{current_pos} < {history['current_max']}. This variant lies outside your defined target area!")

                    g_line, v_line, history = v.get_newlines(line_elements, history, fasta, args.min_depth, args.pl)
                    outfile.write(g_line + '\n')
                    outfile.write(v_line + '\n')

                elif history['current'] != line_elements[0]:
                    # If it isn't,
                    #   write the previous+1 variant position with an end to the max of BED
                    #   Remove previous chrom from acceptable values
                    #   Get new chrom & min value
                    g_line, v_line, history = v.get_newlines(line_elements, history, fasta, args.min_depth, args.pl)
                    # v_line is not used on purpose
                    outfile.write(g_line + '\n')
                    # Remove the BED Regions from the global pool for that chromosome
                    bed_regions = br.filter_chrom_data(line_elements[0], bed_regions)
    # If eof, make sure there are no remaining bed regions
    final_lines = v.get_final_line(history, fasta, args.min_depth, args.pl)
    if final_lines:
        for g_line in final_lines:
            outfile.write(g_line + '\n')
    bed_regions = br.filter_chrom_data(line_elements[0], bed_regions)
    # If any remaining BED regions are left, make sure to print those g_lines as well
    history['current_regions'] = bed_regions
    history['current_max'] = 1
    history['current_min'] = 1
    history['current'] = None
    final_lines = v.get_final_line(history, fasta, args.min_depth, args.pl)
    if final_lines:
        for g_line in final_lines:
            outfile.write(g_line + '\n')
    outfile.close()
