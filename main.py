import pysam
import os
import modules.vcf_handler
import modules.ref_handler
import modules.bam_handler_mpileup
import modules.pileup_creator

def get_pileup_stats(ref, pileup):
    total = 0
    match = 0
    mismatch = 0
    for base in pileup:
        if base == ref:
            match += 1
        else:
            mismatch += 1
        total += 1
    return (total, match, mismatch)


def handle_directory(dir):
    if dir[-1] != '/':
        dir+='/'
    if not os.path.exists(dir):
        os.mkdir(dir)
    return dir


def get_label(genotype_type):
    if genotype_type == "Het":
        return 0
    if genotype_type == "Hom_alt":
        return 1

if '__main__':
    contig = "chr3"
    # site = ""
    site = ":100000-200000"
    vcf_handler = modules.vcf_handler.VCFFileProcessor(
        "/Users/kishwar/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3_split.vcf.gz")
    vcf_handler.populate_dictionary(contig, site)

    ref_handler = modules.ref_handler.RefFileProcessor(
        "/Users/kishwar/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3.fa")
    bam_handler = modules.bam_handler_mpileup.BamProcessor(
        "/Users/kishwar/Kishwar/Whole_chr3_data/illumina/vcf_whole_chr/chr3.bam")

    save_output_dir = handle_directory("/Users/kishwar/Kishwar/pileup-output/Test")
    smry = open(save_output_dir + "summary" + '_' + contig + site.replace(':', '_').replace('-', '_') + ".csv", 'w')

    dict = vcf_handler.get_variant_dictionary()
    for pos in dict.keys():
        for rec in dict[pos]:
            if rec.genotype_class == 'SNP':
                # print(rec)
                pileupcolumns = bam_handler.get_pileupcolumns_aligned_to_a_site(contig, pos-1)
                pileup_object = modules.pileup_creator.PileupProcessor(ref_handler, pileupcolumns, contig, pos-1,
                                                                   rec.type, rec.alt)
                # pileup_object.create_text_pileup(pos-1)
                # pileup_object.create_text_pileup(pos-1)
                image_array, arrayShape = pileup_object.create_image_test(pos-1, image_height=200, image_width=300,
                                                                     ref_band=5, alt=rec.alt)
                file_name = contig + "_" + str(rec.pos)
                pileup_object.save_image_as_png(image_array, save_output_dir, file_name)
                # misc.imsave(self.outputFilename + ".png", pileupArray2d, format="PNG")
                label = get_label(rec.type)
                smry.write(os.path.abspath(save_output_dir + file_name) + ".png," + str(label) + ',' + ','.join(
                    map(str, arrayShape)) + '\n')



    # print("POS\tQUAL\tREF\tALT\tTYPE\tCLASS\tREF\tCOV\tMATCH\tMISM\tPILEUP")
    # dict = vcf_handler.get_variant_dictionary()
    # for pos in dict.keys():
    #     for rec in dict[pos]:
    #         if rec.genotype_class == 'SNP':
    #             pileup_str = bam_handler.get_pileup_of_a_site(contig, pos-1)
    #             # print("HERE", pileup_str)
    #             ref = ref_handler.get_ref_of_region(contig, ":"+str(pos)+"-"+str(pos))
    #             # if ref != rec.ref:
    #                 # raise ValueError("REF DID NOT MATCH FETCHED REF", rec, ref)
    #             (t, m, mm) = get_pileup_stats(ref, pileup_str)
    #             m = int(m*100/t)
    #             mm = int(mm*100/t)
    #             print(rec, "\t" + ref + "\t" + str(t) + "\t" + str(m) + "%\t" + str(mm) + "%\t" + pileup_str)
    #
    #             # bam_handler.test_method(contig, pos-1)
    #             exit()