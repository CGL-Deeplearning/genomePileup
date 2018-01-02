import pysam
import sys


class RefFileProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        try:
            self.FastaFile = pysam.FastaFile(self.file_path)
        except:
            raise IOError("REFERENCE FILE ERROR")

    def get_ref_of_region(self, contig, site):
        return self.FastaFile.fetch(region=contig+site).upper()

