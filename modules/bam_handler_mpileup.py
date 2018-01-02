import pysam
import sys


class BamProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        try:
            self.bamFile = pysam.AlignmentFile(self.file_path, "rb")
        except:
            raise IOError("BAM FILE ERROR")

    def get_pileupcolumns_aligned_to_a_site(self, contig, pos):
        pileupcolumns = self.bamFile.pileup(contig, pos, pos+1)
        return pileupcolumns

    def get_pileup_of_a_site(self, contig, pos):
        pileup_str = ""
        for pileupcolumn in self.bamFile.pileup(contig, pos, pos+1, truncate=True):
            if pileupcolumn.pos < pos:
                continue
            elif pileupcolumn.pos >= pos+1:
                break
            pileup_str += str(pileupcolumn.pos) + " "
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.
                    pileup_str += pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupread.is_del:
                    pileup_str += '*'
        return pileup_str

