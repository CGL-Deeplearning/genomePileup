"""
Implemented by: Kishwar Shafin
Date: 02/01/2018
"""


import pysam
import os
import sys
import numpy
import argparse
from multiprocessing import Process

import modules.vcf_handler
import modules.ref_handler
import modules.bam_handler_mpileup
import modules.pileup_creator


def handle_directory(directory_path):
    """
    Create a directory if doesn't exist
    :param directory_path: path to the directory
    :return: desired directory name
    """
    # if directory has no trailing '/' then add it
    if directory_path[-1] != '/':
        directory_path += '/'
    # if directory doesn't exist then create it
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)

    return directory_path


def get_label(genotype_type):
    """
    Get genotype label for a type of genotype
    :param genotype_type: Genotype in string
    :return: Integer label for the type
    """
    if genotype_type == "Hom":
        return 0
    elif genotype_type == "Het":
        return 1
    elif genotype_type == "Hom_alt":
        return 2


def chunk_it(seq, num):
    """
    Chunk a sequence in N equal segments
    :param seq: Sequence of numbers
    :param num: Number of chunks
    :return: chunked start and end positions
    """
    # find average chunk size
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    # until the end of sequence
    while last < len(seq):
        # append the value to a bin
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def get_alts_in_hom_pileup(pileup_str, ref_base):
    """
    Return possible alts in homozygous cases.
    :param pileup_str: Pileup of a homozygous base.
    :param ref_base: Reference base of that position.
    :return:
    """
    alts = {'A':0, 'C':0, 'G':0, 'T':0}
    for base in pileup_str:
        if base != ref_base and base in alts.keys():
            alts[base] += 1

    return max(alts, key=alts.get), alts[max(alts, key=alts.get)]


def get_odds_for_hom(total_hom, total_het, total_homalt):
    """
    This class will return the odds of generating an image each time we see a hom case
    :param total_hom: Total hom cases present
    :param total_het: Total het cases present
    :param total_homalt: Total homalt cases present
    :return:
    """
    probability_of_seeing_hom = total_hom / (total_hom + total_het + total_homalt)
    odds_of_selecting_hom = 1.0 - probability_of_seeing_hom

    return odds_of_selecting_hom


def generate_pileup(contig, site, bam_file, ref_file, vcf_file, output_dir):
    """
    Generate pileup images from a vcf file
    :param contig: Which contig to fetch ("chr3")
    :param site: Which site to fetch (":100000-200000")
    :param bam_file: Path to the bam alignment file
    :param ref_file: Path to the reference file
    :param vcf_file: Path to the vcf file
    :param output_dir: Output directory, where the image will be saved
    :return:
    """
    # create the vcf handler
    vcf_handler = modules.vcf_handler.VCFFileProcessor(vcf_file)
    # generate dictionary of the region
    vcf_handler.populate_dictionary(contig, site, hom_filter=False)

    # create ref and bam files handler
    ref_handler = modules.ref_handler.RefFileProcessor(ref_file)
    bam_handler = modules.bam_handler_mpileup.BamProcessor(bam_file)

    # create a summary file
    smry = open(output_dir + "summary" + '_' + contig + site.replace(':', '_').replace('-', '_') + ".csv", 'w')

    # get the vcf dictionary of that region
    vcf_dict = vcf_handler.get_variant_dictionary()

    # get the odds of selecting a homozygous case
    total_hom, total_het, total_homalt = vcf_handler.get_genotype_counts()
    odds_of_generating_hom_case = get_odds_for_hom(total_hom, total_het, total_homalt)

    # keep count of how many images of each type is generated
    total_generated_hom, total_generated_het, total_generated_hom_alt = 0, 0, 0

    for pos in vcf_dict.keys():
        in_there = False
        del_there = False
        for rec in vcf_dict[pos]:

            # if rec.genotype_class == 'DEL':
            #     # get pileup columns from bam file
            #     pileup_columns = bam_handler.get_pileup_of_a_site(contig, rec.pos-1)
            # else:
            #     continue
            # if rec.genotype_class == 'DEL' and rec.type == 'Hom_alt':
            #     print(rec)
            if rec.genotype_class == 'IN':
                in_there = True
            if rec.genotype_class == 'DEL':
                del_there = True
        if in_there is True and del_there is True:
            print(contig, pos)
            for rec in vcf_dict[pos]:
                print(rec)

    #         if True:
    #
    #             sys.stderr.write(str(rec)+"\n")
    #             print(rec, end='\t')
    #             # get pileup columns from bam file
    #             pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos - 1)
    #             # create the pileup processor object
    #             pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos - 1,
    #                                                                    rec.type, rec.alt)
    #             # create the image
    #             rgb_image, support_count, unsupport_count = pileup_object.create_image_rgb(pos - 1, image_height=299, image_width=299,
    #                                                                        ref_band=5, alt=rec.alt)
    #
    #             print(support_count, unsupport_count)
    #             file_name = contig + "_" + str(rec.pos)
    #             rgb_image.save(output_dir + file_name + '.png')
    #
    #             if support_count == 0:
    #                 sys.stderr.write('SUPPORT COUNT 0 ENCOUNTERED')
    #                 print('Culprit:\n', rec)
    #                 exit()
    #         # if genotype is SNP then generate image
    #         if rec.genotype_class == 'SNP':
    #             alt = '.'
    #             if rec.type == 'Hom':
    #                 pileup_str = bam_handler.get_pileup_of_a_site(contig, rec.pos-1).split(' ')[1]
    #                 ref_at_pos = ref_handler.get_ref_of_region(contig, ":" + str(rec.pos) + "-" + str(rec.pos))
    #                 alt, mismatches = get_alts_in_hom_pileup(pileup_str, ref_at_pos)
    #                 if mismatches == 0:
    #                     continue
    #
    #             if rec.type == 'Hom' and numpy.random.uniform(0, 1) > odds_of_generating_hom_case:
    #                 continue
    #             elif rec.type == 'Hom':
    #                 rec.alt = alt
    #
    #             total_generated_hom += 1 if rec.type == 'Hom' else 0
    #             total_generated_het += 1 if rec.type == 'Het' else 0
    #             total_generated_hom_alt += 1 if rec.type == 'Hom_alt' else 0
    #
    #             # get pileup columns from bam file
    #             pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos-1)
    #             # create the pileup processor object
    #             pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos-1,
    #                                                                rec.type, rec.alt)
    #             # create the image
    #             image_array, array_shape = pileup_object.create_image_test(pos-1, image_height=299, image_width=299,
    #                                                                        ref_band=5, alt=rec.alt)
    #             # file name for the image and save the image
    #             file_name = contig + "_" + str(rec.pos)
    #             pileup_object.save_image_as_png(image_array, output_dir, file_name)
    #
    #             # label of the image and save the image
    #             label = get_label(rec.type)
    #             smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
    #                 map(str, array_shape)) + '\n')
    #
    #             # report progress
    #             if (total_generated_hom_alt+total_generated_hom+total_generated_het) % 100 == 0:
    #                 total = (total_generated_hom_alt+total_generated_hom+total_generated_het)
    #                 sys.stderr.write(str(total) + ' variants processed in region ' + str(contig) + str(site) + "\n")
    #
    # # print some stats
    # exit()
    # sys.stderr.write('IN REGION: ' + str(contig) + ' ' + site + '\n')
    # sys.stderr.write('TOTAL IN RECORDS:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT\t' + '\n')
    # sys.stderr.write(str(total_hom) + '\t' + str(total_het) + '\t' + str(total_homalt) + '\n')
    #
    # sys.stderr.write('TOTAL GENERATED:\n' + 'HOM\t' + 'HET\t' + 'HOM_ALT' + '\n')
    # sys.stderr.write(str(total_generated_hom) + '\t' + str(total_generated_het) + '\t'
    #                  + str(total_generated_hom_alt) + '\n')


def parallel_pileup_generator(contig, site, bam_file, ref_file, vcf_file, output_dir, threads):
    """
    Generate pileup images from a vcf file using multithreading
    :param contig: Which contig to fetch ("chr3")
    :param site: Which site to fetch (":100000-200000")
    :param bam_file: Path to the bam alignment file
    :param ref_file: Path to the reference file
    :param vcf_file: Path to the vcf file
    :param output_dir: Output directory, where the image will be saved
    :param threads: Number of threads to use for generation
    :return:
    """
    all_positions = []

    for rec in pysam.VariantFile(vcf_file).fetch(region=contig+site):
        all_positions.append(rec.pos)

    all_positions = chunk_it(all_positions, threads)
    starts = [all_positions[i][0] for i in range(len(all_positions))]
    ends = [all_positions[i][-1] for i in range(len(all_positions))]
    args = ()

    for i in range(len(starts)):
        args += ((starts[i], ends[i]),)
        site = ":"+str(starts[i])+"-"+str(ends[i])
        p = Process(target=generate_pileup, args=(contig, site, bam_file, ref_file, vcf_file, output_dir))
        p.start()


if __name__ == '__main__':
    """
    Processes arguments and performs tasks to generate the pileup.
    """

    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file with alignments."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing SNPs and SVs."
    )
    parser.add_argument(
        "--contig",
        type=str,
        default="chr3",
        help="Contig to fetch. Ex: chr3"
    )
    parser.add_argument(
        "--site",
        type=str,
        default="",
        help="Site region. Ex: :100000-200000"
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Name of output directory"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--parallel",
        type=bool,
        default=False,
        help="If true, will use multiple threads."
    )
    FLAGS, not_parsed_flags = parser.parse_known_args()
    # make output directory if not already created
    FLAGS.output_dir = handle_directory(FLAGS.output_dir)

    if FLAGS.parallel:
        parallel_pileup_generator(contig=FLAGS.contig,
                                  site=FLAGS.site,
                                  bam_file=FLAGS.bam,
                                  ref_file=FLAGS.ref,
                                  vcf_file=FLAGS.vcf,
                                  output_dir=FLAGS.output_dir,
                                  threads=FLAGS.max_threads)
    else:
        generate_pileup(contig=FLAGS.contig,
                        site=FLAGS.site,
                        bam_file=FLAGS.bam,
                        ref_file=FLAGS.ref,
                        vcf_file=FLAGS.vcf,
                        output_dir=FLAGS.output_dir)
