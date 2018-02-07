import os
import argparse
import multiprocessing
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


def generate_pileup(bam_file, ref_file, bed_file, output_dir):
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
    smry = open(output_dir + "summary" + '_' + bed_file + ".csv", 'w')

    for rec in all_bed_records:
        print(rec)
        contig, pos, end_pos, ref, alt, genotype = rec.rstrip().split('\t')
        pos = int(pos) + 1

        # get pileup columns from bam file
        pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos)
        # create the pileup processor object
        pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos, genotype, alt)

        # create the image
        image_array, array_shape = pileup_object.create_image(pos - 1, image_height=300, image_width=300, ref_band=5,
                                                              alt=alt, ref=ref)
        # file name for the image and save the image
        file_name = contig + "_" + str(pos) + "_" + str(alt) + "_" + str(genotype)
        pileup_object.save_image_as_png(image_array, output_dir, file_name)

        # label of the image and save the image
        label = genotype
        smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
            map(str, array_shape)) + ',' + str(genotype) + '\n')


def generate_pileup_pl(bam_file, ref_file, bed_records, output_dir, thread_no):
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
    # generate dictionary of the region
    all_bed_records = bed_records

    # create ref and bam files handler
    ref_handler = modules.ref_handler.RefFileProcessor(ref_file)
    bam_handler = modules.bam_handler_mpileup.BamProcessor(bam_file)

    # create a summary file
    smry = open(output_dir + "summary" + '_' + str(thread_no) + ".csv", 'w')

    for rec in all_bed_records:
        contig, pos, end_pos, ref, alt, genotype = rec.rstrip().split('\t')
        pos = int(pos) + 1

        # get pileup columns from bam file
        pileup_columns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos)
        # create the pileup processor object
        pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileup_columns, contig, pos, genotype, alt)

        # create the image
        image_array, array_shape = pileup_object.create_image(pos - 1, image_height=300, image_width=300, ref_band=5,
                                                              alt=alt, ref=ref)
        # file name for the image and save the image
        file_name = contig + "_" + str(pos) + "_" + str(alt) + "_" + str(genotype)
        pileup_object.save_image_as_png(image_array, output_dir, file_name)

        # label of the image and save the image
        label = genotype
        smry.write(os.path.abspath(output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
            map(str, array_shape)) + ',' + str(genotype) + '\n')


def parallel_pileup_generator(bam_file, ref_file, bed_file, output_dir, threads):
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
    all_positions = list()
    bed_handler = modules.bed_handler.BedHandler(bed_file)
    segmented_list_len = 10000
    list_chunks = list()

    for rec in bed_handler.all_bed_records:
        all_positions.append(rec)
        if len(all_positions) >= segmented_list_len:
            list_chunks.append(all_positions)
            all_positions = list()

    if len(all_positions) > 0:
        list_chunks.append(all_positions)

    for i in range(len(list_chunks)):
        p = Process(target=generate_pileup_pl, args=(bam_file, ref_file, list_chunks[i], output_dir, str(i)))
        p.start()
        while True:
            if len(multiprocessing.active_children()) < thread:
                break



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
        parallel_pileup_generator(bam_file=FLAGS.bam,
                                  ref_file=FLAGS.ref,
                                  bed_file=FLAGS.bed,
                                  output_dir=FLAGS.output_dir,
                                  threads=FLAGS.max_threads)
    else:
        generate_pileup(bam_file=FLAGS.bam,
                        ref_file=FLAGS.ref,
                        bed_file=FLAGS.bed,
                        output_dir=FLAGS.output_dir)
