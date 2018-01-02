import pysam
import sys
import numpy as np
from PIL import Image
from scipy import misc

MAX_COLOR_VALUE = 254.0
BASE_QUALITY_CAP = 40.0
MAP_QUALITY_CAP = 60.0
MAP_QUALITY_FILTER = 10.0
class imageChannels:
    def __init__(self, pileup_attributes, ref_base):
        self.pileup_base = pileup_attributes[0]
        self.map_qual = pileup_attributes[1]
        self.base_qual = pileup_attributes[2]
        self.is_rev = pileup_attributes[3]
        self.ref_base = ref_base
        self.is_match = True if self.ref_base == self.pileup_base else False

    @staticmethod
    def get_base_color(base):
        if base == 'A':
            return 250.0
        if base == 'C':
            return 100.0
        if base == 'G':
            return 180.0
        if base == 'T':
            return 30.0
        if base == '*' or 'N':
            return 5.0

    @staticmethod
    def get_base_quality_color(base_quality):
        c_q = min(base_quality, BASE_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / BASE_QUALITY_CAP
        return color

    @staticmethod
    def get_map_quality_color(map_quality):
        c_q = min(map_quality, MAP_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / MAP_QUALITY_CAP
        return color

    @staticmethod
    def get_strand_color(is_rev):
        if is_rev is True:
            return 240
        else:
            return 70

    @staticmethod
    def get_match_ref_color(is_match):
        if is_match is True:
            return MAX_COLOR_VALUE * 0.2
        else:
            return MAX_COLOR_VALUE * 1.0

    @staticmethod
    def get_alt_support_color(is_in_support):
        if is_in_support is True:
            return MAX_COLOR_VALUE * 1.0
        else:
            return MAX_COLOR_VALUE * 0.6

    @staticmethod
    def get_empty_test_channels():
        return [0, 0, 0, 0, 0]

    def get_channels_test(self):
        base_color = self.get_base_color(self.pileup_base)
        base_quality_color = imageChannels.get_base_quality_color(self.base_qual)
        map_quality_color = imageChannels.get_map_quality_color(self.map_qual)
        strand_color = imageChannels.get_strand_color(self.is_rev)
        match_color = imageChannels.get_match_ref_color(self.is_match)

        return [base_color, base_quality_color, map_quality_color, strand_color, match_color]

    @staticmethod
    def get_channels_for_ref_test(base):
        base_color = imageChannels.get_base_color(base)
        base_quality_color = imageChannels.get_base_quality_color(60)
        map_quality_color = imageChannels.get_map_quality_color(60)
        strand_color = imageChannels.get_strand_color(is_rev=False)
        get_match_color = imageChannels.get_match_ref_color(is_match=True)

        return [base_color, base_quality_color, map_quality_color, strand_color, get_match_color]

    # RGB image creator
    @staticmethod
    def get_empty_rgb_channels():
        return [0, 0, 0, 255]

    def get_channels_only_rgb(self):
        if self.ref_base == self.pileup_base:
            return [255, 255, 255, 255]
        elif self.pileup_base == 'A':
            return [255, 0, 0, 255]
        elif self.pileup_base == 'C':
            return [255, 255, 0, 255]
        elif self.pileup_base == 'T':
            return [0, 0, 255, 255]
        elif self.pileup_base == 'G':
            return [0, 255, 0, 255]
        else:
            return [255, 0, 255, 255]

    @staticmethod
    def get_channels_for_ref_only_rgb(base):
        if base == 'A':
            return [255, 0, 0, 255]
        elif base == 'C':
            return [255, 255, 0, 255]
        elif base == 'T':
            return [0, 0, 255, 255]
        elif base == 'G':
            return [0, 255, 0, 255]
        else:
            return [255, 0, 255, 255]



class PileupProcessor:
    def __init__(self, ref_object, pileupcolumns, contig, pos, genotype, alt):
        self.ref_object = ref_object
        self.pileupcolumns = pileupcolumns
        self.contig = contig
        self.pos = pos
        self.genotype = genotype
        self.alt = alt
        # [genomic_position] = [max_insert_length]
        self.insert_length_dictionary = {} # used
        # [read_id] = {{genomic_position}->base}
        self.read_dictionary = {} # used
        # [read_id] = {{genomic_position}->insert_bases}
        self.read_insert_dictionary = {} #used
        # List of Read ids in a genomic position
        self.reads_aligned_to_pos = {}
        # genomic_position_1, genomic_position_2...
        self.position_list = [] # used
        self.leftmost_genomic_position = -1
        self.rightmost_genomic_position = -1
        self.genomic_position_projection = {}
        self.reference_base_projection = {}
        self.ref_sequence = ''
        self.process_pileup()
        self.project_genomic_positions()

    def project_genomic_positions(self):
        ref_seq = self.ref_object.get_ref_of_region(self.contig,
                                                    ":"+str(self.leftmost_genomic_position+1)+ "-"
                                                    + str(self.rightmost_genomic_position+1))
        ref_seq_with_insert = ''
        idx = 0
        for i in range(self.leftmost_genomic_position, self.rightmost_genomic_position+1):
            self.genomic_position_projection[i] = idx
            self.reference_base_projection[i] = ref_seq[i-self.leftmost_genomic_position]
            ref_seq_with_insert += ref_seq[i-self.leftmost_genomic_position]
            idx += 1
            if i in self.insert_length_dictionary:
                ref_seq_with_insert += (self.insert_length_dictionary[i] * '*')
                idx += self.insert_length_dictionary[i]
        self.ref_sequence = ref_seq_with_insert
        return idx

    def length_of_region(self):
        length = 0
        for i in range(self.leftmost_genomic_position, self.rightmost_genomic_position):
            length += 1
            if i in self.insert_length_dictionary:
                length += self.insert_length_dictionary[i]
        return length

    def initialize_dictionaries(self, genomic_position, read_id, is_insert):
        if self.leftmost_genomic_position < 0 or genomic_position < self.leftmost_genomic_position:
            self.leftmost_genomic_position = genomic_position
        if self.rightmost_genomic_position < 0 or genomic_position > self.rightmost_genomic_position:
            self.rightmost_genomic_position = genomic_position

        if genomic_position not in self.reads_aligned_to_pos:
            self.reads_aligned_to_pos[genomic_position] = []

        if read_id not in self.read_dictionary:
            self.read_dictionary[read_id] = {}
            self.read_dictionary[read_id][genomic_position] =''

        if is_insert:
            if genomic_position not in self.insert_length_dictionary:
                self.insert_length_dictionary[genomic_position] = 0
            if read_id not in self.read_insert_dictionary:
                self.read_insert_dictionary[read_id] = {}
                self.read_insert_dictionary[read_id][genomic_position] =''

    def save_info_of_a_position(self, genomic_position, read_id, base, base_qual, map_qual, is_rev, is_in):
        self.initialize_dictionaries(genomic_position, read_id, is_in)

        if is_in is False:
            self.read_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, is_rev)
        else:
            self.read_insert_dictionary[read_id][genomic_position] = (base, base_qual, map_qual, is_rev)
            self.insert_length_dictionary[genomic_position] = max(self.insert_length_dictionary[genomic_position],
                                                                  len(base))

    @staticmethod
    def get_attributes_to_save_indel( pileupcolumn, pileupread):
        insert_start = pileupread.query_position + 1
        insert_end = insert_start + pileupread.indel

        return pileupcolumn.pos, \
               pileupread.alignment.query_name, \
               pileupread.alignment.query_sequence[insert_start:insert_end], \
               pileupread.alignment.query_qualities[insert_start:insert_end], \
               pileupread.alignment.mapping_quality, \
               pileupread.alignment.is_reverse

    @staticmethod
    def get_attributes_to_save(pileupcolumn, pileupread):
        if pileupread.is_del:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name,\
                   '*', \
                   0,\
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse
        else:
            return pileupcolumn.pos, \
                   pileupread.alignment.query_name, \
                   pileupread.alignment.query_sequence[pileupread.query_position],  \
                   pileupread.alignment.query_qualities[pileupread.query_position], \
                   pileupread.alignment.mapping_quality, \
                   pileupread.alignment.is_reverse

    @staticmethod
    def save_image_as_png(pileup_array, save_dir, file_name):
        pileupArray2d = pileup_array.reshape((pileup_array.shape[0], -1))
        misc.imsave(save_dir + file_name + ".png", pileupArray2d, format="PNG")


    def process_pileup(self):
        for pileupcolumn in self.pileupcolumns:
            self.position_list.append(pileupcolumn.pos)
            self.reads_aligned_to_pos[pileupcolumn.pos] = []
            for pileupread in pileupcolumn.pileups:
                self.reads_aligned_to_pos[pileupcolumn.pos].append(pileupread.alignment.query_name)

                if pileupread.indel > 0:
                    gen_pos, read_id, base, base_qual, map_qual, is_rev = \
                        self.get_attributes_to_save_indel(pileupcolumn, pileupread)
                    self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, is_in=True)

                gen_pos, read_id, base, base_qual, map_qual, is_rev = \
                    self.get_attributes_to_save(pileupcolumn, pileupread)
                self.save_info_of_a_position(gen_pos, read_id, base, base_qual, map_qual, is_rev, is_in=False)

    def create_text_pileup(self, query_pos):
        left_most_pos = -1
        for read_id in self.reads_aligned_to_pos[query_pos]:
            read_list = []
            aligned_positions = sorted(self.read_dictionary[read_id].keys())
            if left_most_pos < 0:
                left_most_pos = aligned_positions[0]
            left_most_pos = min(left_most_pos, aligned_positions[0])
            inserts_in_between = sum(self.insert_length_dictionary[val] if val in self.insert_length_dictionary.keys() else 0 for val in range(left_most_pos, aligned_positions[0]))
            padding = aligned_positions[0] - left_most_pos + inserts_in_between
            for pad in range(padding):
                read_list.append(' ')
            for pos in aligned_positions:
                read_list.append((self.read_dictionary[read_id][pos][0]))
                if pos in self.insert_length_dictionary.keys() and self.insert_length_dictionary[pos] > 0:
                    this_has_insert = read_id in self.read_insert_dictionary and pos in self.read_insert_dictionary[read_id]
                    inserted_bases = 0
                    if this_has_insert is True:
                        for base in self.read_insert_dictionary[read_id][pos][0]:
                            read_list.append(base)
                            inserted_bases += 1
                    for i in range(inserted_bases, self.insert_length_dictionary[pos]):
                        read_list.append('*')

            print(''.join(read_list))

    def get_row(self, read_id):
        read_list = {}
        read_insert_list = {}
        aligned_positions = sorted(self.read_dictionary[read_id].keys())
        for pos in aligned_positions:
            read_list[pos] = []
            read_list[pos].append(self.read_dictionary[read_id][pos])
            if pos in self.insert_length_dictionary.keys() and self.insert_length_dictionary[pos] > 0:
                read_insert_list[pos] = []
                inserted_bases = 0
                if read_id in self.read_insert_dictionary and pos in self.read_insert_dictionary[read_id]:
                    inserted_bases = len(self.read_insert_dictionary[read_id][pos][0])
                    read_insert_list[pos].append(self.read_insert_dictionary[read_id][pos])

                for i in range(inserted_bases, self.insert_length_dictionary[pos]):
                    read_attribute_tuple = ('*', [BASE_QUALITY_CAP], self.read_dictionary[read_id][pos][2],
                                            self.read_dictionary[read_id][pos][3])
                    read_insert_list[pos].append(read_attribute_tuple)

        return read_list, read_insert_list

    def get_reference_row_rgb(self, image_width):
        image_row = [imageChannels.get_empty_rgb_channels() for i in range(image_width)]
        for i in range(0, min(len(self.ref_sequence), image_width)):
            image_row[i] = imageChannels.get_channels_for_ref_only_rgb(self.ref_sequence[i])
        return image_row

    def create_image_rgb(self, query_pos, image_height, image_width, ref_band):
        whole_image = []
        for i in range(ref_band):
            whole_image.append(self.get_reference_row_rgb(image_width))

        for read_id in self.reads_aligned_to_pos[query_pos]:
            row_list, row_insert_list = self.get_row(read_id)
            image_row = [imageChannels.get_empty_rgb_channels() for i in range(image_width)]

            for position in sorted(row_list):
                imagechannels_object = imageChannels(row_list[position][0], self.reference_base_projection[position])
                if self.genomic_position_projection[position] < image_width:
                    image_row[self.genomic_position_projection[position]] = imagechannels_object.get_channels_only_rgb()

                if position in row_insert_list.keys():
                    insert_ref = 0
                    for bases in row_insert_list[position]:
                        for base_idx in range(len(bases[0])):
                            insert_ref += 1
                            attribute_tuple = (bases[0][base_idx], bases[1][base_idx], bases[2], bases[3])
                            imagechannels_object = imageChannels(attribute_tuple, '*')
                            if self.genomic_position_projection[position] + insert_ref < image_width:
                                image_row[self.genomic_position_projection[position] + insert_ref] = \
                                    imagechannels_object.get_channels_only_rgb()

            whole_image.append(image_row)

        empty_rows_to_add = image_height - len(whole_image)
        for i in range(empty_rows_to_add):
            whole_image.append([imageChannels.get_empty_rgb_channels() for i in range(image_width)])

        image_array = np.array(whole_image).astype(np.uint8)
        img = Image.fromarray(image_array)
        return img

    # TEST FIVE CHANNELS
    def get_reference_row_test(self, image_width):
        image_row = [imageChannels.get_empty_test_channels() for i in range(image_width)]
        for i in range(0, min(len(self.ref_sequence), image_width)):
            image_row[i] = imageChannels.get_channels_for_ref_test(self.ref_sequence[i])
        return image_row

    def create_image_test(self, query_pos, image_height, image_width, ref_band, alt):
        whole_image = []
        for i in range(ref_band):
            whole_image.append(self.get_reference_row_test(image_width))

        for read_id in self.reads_aligned_to_pos[query_pos]:
            row_list, row_insert_list = self.get_row(read_id)
            image_row = [imageChannels.get_empty_test_channels() for i in range(image_width)]
            filter_row = False
            for position in sorted(row_list):
                if row_list[position][0][2] < MAP_QUALITY_FILTER:
                    filter_row = True
                    break
                imagechannels_object = imageChannels(row_list[position][0], self.reference_base_projection[position])
                if self.genomic_position_projection[position] < image_width:
                    image_row[self.genomic_position_projection[position]] = imagechannels_object.get_channels_test()

                if position in row_insert_list.keys():
                    insert_ref = 0
                    for bases in row_insert_list[position]:
                        for base_idx in range(len(bases[0])):
                            insert_ref += 1
                            attribute_tuple = (bases[0][base_idx], bases[1][base_idx], bases[2], bases[3])
                            imagechannels_object = imageChannels(attribute_tuple, '*')
                            if self.genomic_position_projection[position] + insert_ref < image_width:
                                image_row[self.genomic_position_projection[position] + insert_ref] = \
                                    imagechannels_object.get_channels_test()
            if filter_row is False and len(whole_image) < image_height:
                whole_image.append(image_row)

        empty_rows_to_add = image_height - len(whole_image)
        for i in range(empty_rows_to_add):
            whole_image.append([imageChannels.get_empty_test_channels() for i in range(image_width)])

        image_array = np.array(whole_image).astype(np.uint8)
        return image_array, image_array.shape

