"""
Implemented by: Kishwar Shafin
Date: 02/01/2018
"""
import pysam

"""
This class handles fasta file using pysam API. Two methods mostly handles the FASTA interaction:
- __init__: create FastaFile object given file path to a bam file
- get_ref_of_region: return a string containing reference of a site
"""

class RefFileProcessor:
    def __init__(self, file_path):
        """
        create FastaFile object given file path to a bam file
        :param file_path: Path to FASTA file
        """
        self.file_path = file_path
        try:
            self.FastaFile = pysam.FastaFile(self.file_path)
        except:
            raise IOError("REFERENCE FILE ERROR")

    def get_ref_of_region(self, contig, site):
        """
        Return a string containing reference of a site
        :param contig: Contig [ex chr3]
        :param site: Site [ex 100000-200000]
        :return:
        """
        return self.FastaFile.fetch(region=contig+site).upper()

