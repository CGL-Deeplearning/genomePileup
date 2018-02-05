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

import modules.bed_handler
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


def generate_pileup(contig, site, bam_file, ref_file, bed_file, output_dir):
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
    bed_handler = modules.bed_handler.BedHandler(bed_file)

    # generate dictionary of the region
    all_bed_records = bed_handler.all_bed_records

    # create ref and bam files handler
    ref_handler = modules.ref_handler.RefFileProcessor(ref_file)
    bam_handler = modules.bam_handler_mpileup.BamProcessor(bam_file)

    # create a summary file
    smry = open(output_dir + "summary" + '_' + contig + site.replace(':', '_').replace('-', '_') + ".csv", 'w')

    for rec in all_bed_records:
        chr_name, pos, end_pos, ref, alt, genotype = rec.rstrip().split('\t')
        pos = int(pos) + 1

        # get pileup columns from bam file
        pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos)
        # create the pileup processor object
        pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos, genotype, alt)

        # create the image
        image_array, array_shape = pileup_object.create_image(pos - 1, image_height=300, image_width=300, ref_band=5,
                                                              alt=alt, ref=ref)
        # file name for the image and save the image
        file_name = contig + "_" + str(pos)
        pileup_object.save_image_as_png(image_array, output_dir, file_name)

        # label of the image and save the image
        label = genotype
        smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
            map(str, array_shape)) + ',' + str(genotype) + '\n')


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
        "--bed",
        type=str,
        required=True,
        help="Bed file containing labeled records."
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
                        bed_file=FLAGS.bed,
                        output_dir=FLAGS.output_dir)