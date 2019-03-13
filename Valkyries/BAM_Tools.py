"""
This is a collection of functions to make manipulation of BAM files simpler.

BAM_Tools v0.6.0
    March 16, 2018
    Dennis A. Simpson
    Gave the user the ability to control the compression level of the sorted BAM files.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC
@copyright: 2018
"""

import ntpath
from collections import defaultdict, OrderedDict
from copy import deepcopy
import subprocess
import pathlib
import itertools
import pathos
import pybedtools
import pysam
from natsort import natsort
import Valkyries.Utilities as Utilities
import Valkyries.Tool_Box as Tool_Box

__author__ = "Dennis A. Simpson"
__version__ = "0.6.1"


def _byte_array_to_string(sequence):
    if isinstance(sequence, str):
        return sequence
    else:
        return str(sequence.decode("utf-8"))


class AlignWriter:
    class _NullWriter(object):
        def write(self, family, paired_align, alignment):
            pass

        def close(self, log=None):
            pass

    NULL = _NullWriter()

    def __init__(self, header, bam_path, tags=None):
        if tags is None:
            self._tags = []
        else:
            self._tags = sorted(tags)

        new_header = deepcopy(header)
        if 'CO' not in new_header:
            new_header['CO'] = []

        new_header['CO'].extend([tag.header_comment for tag in self._tags])
        self._bam_path = bam_path
        self._bam_file = pysam.AlignmentFile(bam_path, "wb", header=new_header)

    @property
    def bam_file_path(self):
        return self._bam_path

    def _add_bam_tags(self, family, paired_align, alignment_object):
        for tag in self._tags:
            tag.set_tag(family, paired_align, alignment_object)

    def write(self, family, paired_align, alignment_object):
        self._add_bam_tags(family, paired_align, alignment_object)
        self._bam_file.write(alignment_object.pysam_align_segment)

    def close(self, log=None):
        self._bam_file.close()


class LoggingWriter:

    class UnplacedFamily:
        def __init__(self):
            self.filter_value = 'unplaced'
            self.umi_sequence = -1

    UNPLACED_FAMILY = UnplacedFamily()

    def __init__(self, base_writer, log, chromosome):
        self._base_writer = base_writer
        self._log = log
        self._align_filter_stats = defaultdict(int)
        self._family_filter_stats = defaultdict(set)
        self.chromosome = chromosome

    def write(self, family, paired_align, alignment):
        if not family:
            family = LoggingWriter.UNPLACED_FAMILY

        self._align_filter_stats[(family.filter_value, alignment.filter_value)] += 1
        self._family_filter_stats[family.filter_value].add(family.umi_sequence)

        if family == LoggingWriter.UNPLACED_FAMILY:
            self._base_writer.write(None, paired_align, alignment)
        else:
            self._base_writer.write(family, paired_align, alignment)

    @staticmethod
    def _log_line(text, count, total, filter_name):
        line = '{:.2f}% ({}/{}) {}: {}'
        return line.format(100 * count / total, count, total, text, filter_name)

    @property
    def _unplaced_aligns(self):
        unplaced_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value and align_filter:
                unplaced_aligns[align_filter] = cnt

        return OrderedDict(Tool_Box.sort_dict(unplaced_aligns))

    @staticmethod
    def _discarded_filter_value(fam_filter, align_filter):
        if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value:
            return None
        else:
            filter_values = []
            if fam_filter:
                filter_values.append(fam_filter)
            if align_filter:
                filter_values.append(align_filter)
            return "; ".join(filter_values)

    @property
    def _discarded_aligns(self):
        discarded_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            filter_value = LoggingWriter._discarded_filter_value(fam_filter,
                                                                 align_filter)
            if filter_value:
                discarded_aligns[filter_value] = cnt

        return OrderedDict(Tool_Box.sort_dict(discarded_aligns))

    @property
    def _family_stats(self):
        family_filter_stats = dict(self._family_filter_stats)
        family_filter_stats.pop(LoggingWriter.UNPLACED_FAMILY.filter_value, None)
        included_count = len(family_filter_stats.pop(None, []))
        discarded_count = 0
        filter_counts = OrderedDict()

        for name, fam_ids in family_filter_stats.items():
            align_count = len(fam_ids)
            discarded_count += align_count
            filter_counts[name] = align_count

        total_count = included_count + discarded_count

        return included_count, total_count, OrderedDict(Tool_Box.sort_dict(filter_counts))

    @property
    def _align_stats(self):
        included_filter = (None, None)
        included_count = self._align_filter_stats[included_filter]
        excluded_count = sum([count for fam_align_filter,
                              count in self._align_filter_stats.items() if fam_align_filter != included_filter])
        total_count = included_count + excluded_count
        return included_count, excluded_count, total_count

    @staticmethod
    def _percent_stat_str(count, total):
        return '{:.2f}% ({}/{})'.format(100 * count / total, count, total)

    @staticmethod
    def _log_filter_counts(filter_counts, log_method, msg_format, total, chromosome):
        for name, count in filter_counts.items():
            percent = LoggingWriter._percent_stat_str(count, total)
            log_method(msg_format.format(chromosome=chromosome, filter_name=name, percent_stat=percent))

    def _log_results(self):
        (included_align_count, excluded_align_count, tot_align_count) = self._align_stats
        (included_fam_count, total_fam_count, discarded_fam_filter_counts) = self._family_stats

        self._log.info('{}: {} alignments unplaced or discarded',
                       self.chromosome, LoggingWriter._percent_stat_str(excluded_align_count, tot_align_count))

        LoggingWriter._log_filter_counts(self._unplaced_aligns, self._log.debug,
                                         '{chromosome}: alignments unplaced: {percent_stat} {filter_name}',
                                         tot_align_count, self.chromosome)

        LoggingWriter._log_filter_counts(self._discarded_aligns, self._log.debug,
                                         '{chromosome} alignments discarded: {percent_stat} {filter_name}',
                                         tot_align_count, self.chromosome)

        LoggingWriter._log_filter_counts(discarded_fam_filter_counts, self._log.info,
                                         '{chromosome} families discarded: {percent_stat} {filter_name}',
                                         total_fam_count, chromosome=self.chromosome)

        percent_stat = LoggingWriter._percent_stat_str(included_align_count, tot_align_count)
        self._log.info('{}: {} alignments included in {} families', self.chromosome, percent_stat, included_fam_count)

        if included_align_count == 0:
            self._log.warning("{}: No alignments passed filters. (Was input BAM downsampled?)".format(self.chromosome))
        else:
            percent_dedup = 100 * (1 - (included_fam_count / included_align_count))
            msg = '{} {:.2f}% deduplication rate (1 - {} families/{} included alignments)'

            self._log.info(msg, self.chromosome, percent_dedup, included_fam_count, included_align_count)

    def close(self, log=None):
        if self._align_filter_stats:
            self._log_results()
        self._base_writer.close(log)


class BamTag:

    class _NullObject(object):
        """Returns None for all method calls"""

        def __init__(self):
            self.included_pair_count = None
            self.filter_value = None
            self.umi_sequence = None
            self.umt = lambda *args: None
            self.is_consensus_template = lambda *args: None
            self.positions = lambda *args: None
            self.cigars = lambda *args: None

    _NULL_OBJECT = _NullObject()

    def __init__(self, tag_name, tag_type, description, get_value, analysis):
        self._tag_name = tag_name
        self._tag_type = tag_type
        self._get_value = get_value
        self._description = description

        self.header_comment = "{}\tBAM tag\t{}: {}".format(analysis, tag_name, description)

    def __lt__(self, other):
        return (self._tag_name, self._description) < (other._tag_name, other._description)

    def set_tag(self, family, paired_align, align_obj):
        family = family if family else BamTag._NULL_OBJECT
        paired_align = paired_align if paired_align else BamTag._NULL_OBJECT
        value = self._get_value(family, paired_align, align_obj)
        align_obj.set_tag(self._tag_name, value, self._tag_type)


class BamFlag:
    PAIRED = 1
    PROPER_PAIR = 2
    UNMAP = 4
    MUNMAP = 8
    REVERSE = 16
    MREVERSE = 32
    READ1 = 64
    READ2 = 128
    SECONDARY = 256
    QCFAIL = 512
    DUP = 1024
    SUPPLEMENTARY = 2048


class PairedAlignment:
    """Represents the left and right align pairs from an single sequence."""
    def __init__(self, left_alignment, right_alignment):

        if left_alignment.query_name.count(":") > 7:
            left_query_name = ":".join(left_alignment.query_name.split(":")[:-4])
            right_query_name = ":".join(right_alignment.query_name.split(":")[:-4])
        else:
            left_query_name = left_alignment.query_name
            right_query_name = right_alignment.query_name

        if left_query_name != right_query_name:
            msg = 'Inconsistent query names ({} != {})'
            raise ValueError(msg.format(left_query_name, right_query_name))

        self.query_name = left_alignment.query_name
        self.left = left_alignment
        self.right = right_alignment

        left_umt = left_alignment.query_name.split(":")[0].split("|")[1]
        right_umt = right_alignment.query_name.split(":")[0].split("|")[1]

        if len(right_alignment.query_name.split(":")[0].split("|")) == 3:
            right_umt = right_alignment.query_name.split(":")[0].split("|")[2]
            left_umt = left_alignment.query_name.split(":")[0].split("|")[1]

        self.umt = (left_umt, right_umt)
        self._tag_length = len(left_umt)

    @property
    def filter_value(self):
        if self.left.filter_value or self.right.filter_value:
            return self.left.filter_value, self.right.filter_value
        else:
            return None

    def cigars(self, format_string=None):
        """
        This is the CIGAR used in the tags.
        :param format_string:
        :return:
        """
        if format_string:
            return format_string.format(left=self.left.cigarstring, right=self.right.cigarstring)

        else:
            return self.left.cigarstring, self.right.cigarstring

    def positions(self, format_string=None):
        """
        Make Samtools positions.
        :param format_string:
        :return:
        """
        left_value = self.left.reference_start + 1
        right_value = self.right.reference_end + 1

        if format_string:
            return format_string.format(left=left_value, right=right_value)
        else:
            return left_value, right_value

    # def replace_umt(self, umt):
    #     if not (umt[0] or umt[1]) or (len(umt[0]) != self._tag_length) or (len(umt[1]) != self._tag_length):
    #         msg = "Each UMT must match tag_length ({})"
    #         raise ValueError(msg.format(self._tag_length))

        # left_qual = self.left.query_qualities
        # right_qual = self.right.query_qualities
        #
        # left_query_frag = self.left.query_sequence[len(umt[0]):]
        # left_query_frag_str = _byte_array_to_string(left_query_frag)
        #
        # self.left.query_sequence = umt[0] + left_query_frag_str
        # right_query_frag = self.right.query_sequence[:-len(umt[1])]
        # right_query_frag_str = _byte_array_to_string(right_query_frag)
        # self.right.query_sequence = right_query_frag_str + umt[1]
        # self.umt = umt
        # self.left.query_qualities = left_qual
        # self.right.query_qualities = right_qual

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(self.left) * hash(self.right)

    def __repr__(self):
        return ("Pair({}|{}|{}, "
                "{}|{}|{})").format(self.left.query_name,
                                    self.left.reference_start,
                                    self.left.query_sequence,
                                    self.right.query_name,
                                    self.right.reference_start,
                                    self.right.query_sequence)


class MimirAlign:
    # cgates: FYI, you can use dynamic delegation via __setattr__ and
    # __getattr__ but it's awkward and about twice as slow
    def __init__(self, pysam_align_segment, filter_value=None):
        self.pysam_align_segment = pysam_align_segment
        self.filter_value = filter_value

    def __eq__(self, other):
        return other.__dict__ == self.__dict__

    # cgates: the native pysam hashing is not performant for ultradeep pileups
    def __hash__(self):
        return hash(self.filter_value) ^ \
               hash(self.pysam_align_segment.query_name) ^ \
               self.pysam_align_segment.reference_start

    @property
    def cigarstring(self):
        return self.pysam_align_segment.cigarstring

    @cigarstring.setter
    def cigarstring(self, value):
        self.pysam_align_segment.cigarstring = value

    @property
    def flag(self):
        return self.pysam_align_segment.flag

    @flag.setter
    def flag(self, value):
        self.pysam_align_segment.flag = value

    def get_tag(self, name, with_value_type=False):
        return self.pysam_align_segment.get_tag(name, with_value_type)

    def get_tags(self, with_value_type=False):
        return self.pysam_align_segment.get_tags(with_value_type)

    @property
    def mapping_quality(self):
        return self.pysam_align_segment.mapping_quality

    @mapping_quality.setter
    def mapping_quality(self, value):
        self.pysam_align_segment.mapping_quality = value

    @property
    def next_reference_start(self):
        return self.pysam_align_segment.next_reference_start

    @next_reference_start.setter
    def next_reference_start(self, value):
        self.pysam_align_segment.next_reference_start = value

    @property
    def orientation(self):
        if self.reference_start < self.next_reference_start:
            return 'left'
        elif self.reference_start > self.next_reference_start:
            return 'right'
        else:
            return 'neither'

    @property
    def query_name(self):
        return self.pysam_align_segment.query_name

    @query_name.setter
    def query_name(self, value):
        self.pysam_align_segment.query_name = value

    @property
    def query_sequence(self):
        return self.pysam_align_segment.query_sequence

    @query_sequence.setter
    def query_sequence(self, value):
        self.pysam_align_segment.query_sequence = value

    @property
    def query_qualities(self):
        return self.pysam_align_segment.query_qualities

    @query_qualities.setter
    def query_qualities(self, value):
        self.pysam_align_segment.query_qualities = value

    @property
    def reference_end(self):
        return self.pysam_align_segment.reference_end

    @property
    def reference_id(self):
        return self.pysam_align_segment.reference_id

    @reference_id.setter
    def reference_id(self, value):
        self.pysam_align_segment.reference_id = value

    @property
    def reference_name(self):
        return self.pysam_align_segment.reference_name

    @property
    def reference_start(self):
        return self.pysam_align_segment.reference_start

    @reference_start.setter
    def reference_start(self, value):
        self.pysam_align_segment.reference_start = value

    def set_tag(self, tag_name, tag_value, value_type):
        self.pysam_align_segment.set_tag(tag_name, tag_value, value_type)

    @property
    def template_length(self):
        return self.pysam_align_segment.template_length

    @template_length.setter
    def template_length(self, value):
        self.pysam_align_segment.template_length = value


def filter_alignments(pysam_alignments, excluded_writer=AlignWriter.NULL):
    filters = {'cigar unavailable':
                    lambda a: a.cigarstring is None,
               'mapping quality < 1':
                    lambda a: a.mapping_quality < 1,
               'not in proper pair':
                    lambda a: a.flag & BamFlag.PROPER_PAIR == 0,
               'qc failed':
                    lambda a: a.flag & BamFlag.QCFAIL != 0,
               'secondary alignment':
                    lambda a: a.flag & BamFlag.SECONDARY != 0,
               'supplementary alignment':
                    lambda a: a.flag & BamFlag.SUPPLEMENTARY != 0,
               }

    generator = Utilities.FilteredGenerator(filters)
    for pysam_align, filter_value in generator.filter(pysam_alignments):
        alignment = MimirAlign(pysam_align, filter_value)

        if filter_value:
            excluded_writer.write(family=None, paired_align=None, alignment=alignment)
        else:
            yield alignment


def alignment_file(filename, mode, template=None):
    return pysam.AlignmentFile(filename, mode, template)


def sort(input_bam_filepath, threads, compression_level):

    sorted_bamfile = input_bam_filepath.replace(".bam", "_sorted.bam")
    pysam.samtools.sort("-@", str(threads), "-l", compression_level, '-o', sorted_bamfile, input_bam_filepath,
                        catch_stdout=False)

    Tool_Box.delete([input_bam_filepath])
    return sorted_bamfile


def index(bam_filepath):
    pysam.samtools.index(bam_filepath, catch_stdout=False)


def merge_bam(args, log, file_list, merged_bamfile=None):
    """
    Merge BAM files that resulted from alignment of the temporary FASTQ files.  Then sort and index merged BAM.
    :return:
    """

    log.info("Merging tmp BAM files.")
    if not merged_bamfile:
        merged_bamfile = "{0}{1}_merged.bam".format(args.Working_Folder, args.Job_Name)
    bam_files = " ".join(file_list)

    log.debug("Cannot figure out how to get pysam.merge() to work so I use Samtools directly.")
    cmd = "samtools merge -f -O BAM {0} {1}".format(merged_bamfile, bam_files)
    subprocess.run([cmd], shell=True)
    # pysam.samtools.merge("-f", merged_bamfile, bam_files, catch_stdout=False)

    log.info("Temporary BAM files merged.  Sort and index merged BAM file.")
    sorted_bamfile = sort(merged_bamfile, int(args.Spawn), args.Compression_Level)
    index(sorted_bamfile)
    log.info("BAM file sorted and indexed.  Deleting Temporary BAM files..")
    Tool_Box.delete(file_list)

    return sorted_bamfile


def build_writer(bamfile, output_bam, tags, analysis, version):
    if not output_bam:
        return AlignWriter.NULL
    else:
        input_bam = alignment_file(bamfile, 'rb')
        hd = input_bam.header
        input_bam.close()
        header = hd.to_dict()

        HEADER_PG_KEY = 'PG'
        pg_headers = header.get(HEADER_PG_KEY, [])
        pg_headers.append({"ID": analysis, "PN": analysis, "VN": version})
        header[HEADER_PG_KEY] = pg_headers

        return AlignWriter(header, output_bam, tags)


def total_align_count(input_bam, chromosome=None):
    """Returns count of all mapped alignments in input BAM (based on index)"""
    count = 0

    for line in pysam.samtools.idxstats(input_bam).split('\n'):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')

            if chrom != '*' and not chromosome:
                count += int(mapped) + int(unmapped)
            elif chrom == chromosome:
                count += int(mapped) + int(unmapped)
    return count


def build_bam_tags(analysis):

    def combine_filters(fam, pair, align):
        filters = [x.filter_value for x in [fam, align] if x and x.filter_value]

        if filters:
            return ";".join(filters).replace('; ', ';')
        else:
            return None

    boolean_tag_value = {True: 1}
    tags = [
        BamTag("X0", "Z", "filter (why the alignment was excluded)", combine_filters, analysis),
        BamTag("X1", "Z", "leftmost~rightmost matched pair positions",
               lambda fam, pair, align: pair.positions('{left}~{right}'), analysis),
        BamTag("X2", "Z", "L~R CIGARs",
               lambda fam, pair, align: pair.cigars('{left}~{right}'), analysis),
        BamTag("X3", "i", "unique identifier for this alignment family",
               lambda fam, pair, align: fam.umi_sequence, analysis),
        BamTag("X4", "Z", "L~R UMT barcodes for this alignment family; because of fuzzy matching the family UMT may "
                          "be distinct from the UMT of the original alignment",
               lambda fam, pair, align: fam.umt('{left}~{right}'), analysis),
        BamTag("X5", "i", "family size (number of align pairs in this family)",
               lambda fam, pair, align: fam.included_pair_count, analysis),
        BamTag("X6", "i", "presence of this tag signals that this alignment would be the template for the consensus "
                          "alignment",
               lambda fam, pair, align: boolean_tag_value.get(fam.is_consensus_template(align), None), analysis)]
    return tags


class BamTools:
    def __init__(self, args, log, index_list, file_info=None):
        self.log = log
        self.index_list = index_list
        self.file_info = file_info
        self.args = args
        self.bamfile = None
        self.bamfile_list = None
        self.coverage_dict = {"chrom": {"cell": "tuple"}}

    def build_dictionary(self, input_file):
        """
        Get the number of reads aligned for each chromosome and put this into a dictionary.  The "*" in position 0 are
        the total unmapped reads.  The idxstats format is region (Chr usually), region length, mapped, unmapped.
        :return:
        """

        if len(input_file) < 3:
            self.log.error("BAM file parameter missing from options file.")
            raise SystemExit(1)

        elif not pathlib.Path(input_file).is_file():
            self.log.error("BAM file {0} not found".format(input_file))
            raise SystemExit(1)

        algn_dict = {}

        for item in pysam.samtools.idxstats(input_file, split_lines=True):
            if not item.split("\t")[0] == "*":
                algn_dict[item.split("\t")[0]] = (int(item.split("\t")[1]), int(item.split("\t")[2]))

        return algn_dict

    def cell_coverage(self):

        cell_label_string = ["Chrom\tCoverage"]
        cell_breadth_string = []
        cell_depth_string = []
        cell_idxstats_string = []
        cell_idxstats_dict = Tool_Box.VivifiedDictionary()
        first_line = True
        p = pathos.multiprocessing.Pool(int(self.args.Spawn))

        for bamfile in natsort.natsorted(self.bamfile_list):
            cell_name = ntpath.basename(bamfile)
            cell_algn_list = list(self.build_dictionary(bamfile).items())
            cover_tool = Tool_Box.CoverageCalculator(bamfile)

            tmp_dict_list = p.starmap(self.__sub_cell_coverage, zip(itertools.repeat(cover_tool), cell_algn_list))

            for d in tmp_dict_list:
                self.coverage_dict[cell_name].update(d)
            del tmp_dict_list

            cell_idxstats_dict[cell_name] = self.build_dictionary(bamfile)

        for cell_name in natsort.natsorted(self.coverage_dict):
            cell_label_string += "\t{0}".format(cell_name)

            for chrom in natsort.natsorted(self.coverage_dict[cell_name]):
                if first_line:
                    cell_breadth_string += "{0}\tBreadth".format(chrom)
                    cell_depth_string += "{0}\tDepth".format(chrom)
                    cell_idxstats_string += "{0}\tMapped_Reads".format(chrom)
                    first_line = False
                else:
                    cell_breadth_string += "\t{0}".format(self.coverage_dict[cell_name][chrom][0])
                    cell_depth_string += "\t{0}".format(self.coverage_dict[cell_name][chrom][1])
                    cell_idxstats_string += "\t{0}".format(cell_idxstats_dict[cell_name][chrom][0])

            first_line = False
            cell_breadth_string += "\n"
            cell_depth_string += "\n"
            cell_idxstats_string += "\n"

        cell_label_string += "\n"

        outfile = open("{0}{1}_cell_coverage.txt".format(self.args.Working_Folder, self.args.Job_Name), "w")
        outfile.write("{0}{1}{2}".format(cell_label_string, cell_breadth_string, cell_depth_string))
        outfile.close()

    @staticmethod
    def __sub_cell_coverage(coverage_tool, region):
        coverage_dict = {}

        coverage_depth, coverage_breadth = coverage_tool.coverage(region)
        coverage_dict[region[0]] = {coverage_tool.cell_name: {coverage_breadth, coverage_depth}}

        return coverage_dict

    def merge_bam(self):
        """
        Merge BAM files that resulted from alignment of the temporary FASTQ files.  Fills the self.bamfile name.
        :return:
        """

        self.log.info("Merging tmp BAM files.")
        file_info = self.file_info

        self.bamfile = "{0}{1}_merged.bam".format(self.args.Working_Folder, self.args.Job_Name)

        bam_files = " ".join(self.file_info.bam_file_list)
        cmd = "samtools merge {0} {1}".format(self.bamfile, bam_files)
        subprocess.run([cmd], shell=True)

        self.log.info("Temporary BAM files merged.  Sort and index merged BAM file.")
        index(sort(self.bamfile, int(self.args.Spawn), self.args.Compression_Level))
        self.log.info("BAM file sorted and indexed.  Deleting Temporary BAM files..")

        # delete temporary BAM files used in the merge.
        Tool_Box.delete(file_info.bam_file_list)

    def demultiplex(self):
        """
        Take a multiplexed BAM file and create demultiplexed BAM files based on the index list.
        :return:
        """

        self.log.info("-->Demultiplexing BAM File")
        demultiplex_dict = {}
        barcode_list = []
        tmp_bamfile_list = []

        # Grab a barcode from index list so we can get the length.
        for item in self.index_list:
            d5 = Sequence_Magic.rcomp(item[0].split("+")[0][1:])
            d7 = item[0].split("+")[1][:-1]
            barcode = "{0}+{1}".format(d7, d5)
            barcode_length = len(barcode)
            break

        index_list.append(("Unidentified", "Unidentified"))

        # Initialize the output file.
        outfile = "{0}{1}_BAM_QualityAssessment.txt".format(self.args.Working_Folder, self.args.Job_Name)
        quality_outfile = open(outfile, "w")

        # Build the header string and write it to the output file.
        print("\tBuilding Header String for {0}{1}_BAM_QualityAssessment.txt"
              .format(self.args.Working_Folder, self.args.Job_Name))

        write_string1 = "Barcode\tName\tReads"

        for i in range(int(self.args.Mismatch) + 1):
            write_string1 += "\t{0}_Mismatch".format(i)

        write_string1 += "\tN_Count"

        for i in range(barcode_length):
            write_string1 += "\tPosition_{0}".format(i + 1)

        write_string1 += "\n"

        quality_outfile.write(write_string1)

        print("\tHeaders built and file writen.")

        # Begin processing BAM File.
        multiplexed_bamfile = pysam.AlignmentFile(self.bamfile, 'rb')

        print("\t-->Create List of BAM file names and Dictionary of BAM files.")
        for item in index_list:
            if item[0] == "Unidentified":
                barcode = "Unidentified"
            else:
                d5 = Sequence_Magic.rcomp(item[0].split("+")[0][1:])
                d7 = item[0].split("+")[1][:-1]
                barcode = "{0}+{1}".format(d7, d5)

            barcode_list.append((barcode, item[1]))

            bamfile_name = "{0}{1}_{2}.bam".format(self.file_info.fastq1.filepath, self.args.Job_Name, item[1])
            tmp_bamfile_list.append(bamfile_name)
            demultiplex_dict[barcode] = (pysam.AlignmentFile(bamfile_name, 'wb', template=multiplexed_bamfile),
                                         [0] * (int(self.args.Mismatch) + barcode_length + 3))

        print("-->BAM file list and file dictionary initialized.  Begin iterating BAM file.")

        read_count = 0
        for read in multiplexed_bamfile.fetch(until_eof=True):
            if read_count % int(self.args.prog_check) == 0:
                print("\tProcessed {0} segments in BAM File {1}.".format(read_count, multiplexed_bamfile))

            mismatch_found = False
            read_count += 1
            offset = 1
            read_barcode = read.qname.split("|")[1].split(":")[3]

            print(read_barcode)
            for index_query in demultiplex_dict:

                mismatch1 = Sequence_Magic.match_maker(index_query, read_barcode, self.args.Mismatch)

                if mismatch1:
                    mismatch_found = True
                    demultiplex_dict[index_query][0].write(read)

                    # Total reads for barcode
                    demultiplex_dict[index_query][1][0] += 1

                    for i in range(int(self.args.Mismatch) + 1):
                        if i == mismatch1[1]:
                            demultiplex_dict[index_query][1][offset] += 1
                        offset += 1

                    if read_barcode.count("N") > 0:
                        demultiplex_dict[index_query][1][offset] += 1

                    offset += 1

                    for i in range(barcode_length):
                        if read_barcode[i] == "N":
                            demultiplex_dict[index_query][1][i + offset] += 1

                    break

            if not mismatch_found:
                offset = int(self.args.Mismatch) + 2
                demultiplex_dict["Unidentified"][0].write(read)
                demultiplex_dict["Unidentified"][1][0] += 1

                if read_barcode.count("N") > 0:
                    demultiplex_dict["Unidentified"][1][offset] += 1

                offset += 1

                for i in range(barcode_length):
                    if read_barcode[i] == "N":
                        demultiplex_dict["Unidentified"][1][i + offset] += 1

        for key, name in barcode_list:
            demultiplex_dict[key][0].close()
            quality_outfile.write("{0}\t{1}\t{2}\n".format(key, name, '\t'.join(map(str, demultiplex_dict[key][1]))))

        quality_outfile.close()
        multiplexed_bamfile.close()
        demultiplex_dict.clear()

        # Sort and index the demultiplexed BAM files.
        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        self.bamfile_list = \
            (p.starmap(self.bam_sort, zip(tmp_bamfile_list,
                                          itertools.repeat(int(self.args.Spawn)),
                                          itertools.repeat(int(self.args.Compression_Level)))))

        return self.bamfile_list

    def bam_to_bed(self):
        """
        Calls system program to convert BAM file(s) to BED file(s).
        """
        if self.bamfile_list is None:
            raise SystemExit("\033[1;31m***BAM File List is Empty***\033[m")

        count = 0
        bedfile_list = []
        print("Begin converting {0} BAM files to BED files".format(len(self.bamfile_list)))

        for bamfile in self.bamfile_list:
            if not pathlib.Path(bamfile).is_file():
                print("\t\033[1;31mCaution:\033[m File {0} not found.  Results in question".format(bamfile))
            else:
                if pathlib.Path(bamfile).stat().st_size > 500:  # Exclude empty BAM files.
                    print("   -->Converting file \u001b[46;1m{0}\033[m of \u001b[46;1m{1}\033[m.\n\t\t{2}"
                          .format(count + 1, len(self.bamfile_list), bamfile))

                    Tool_Box.debug_messenger("Testing pybedtools for bed_to_bam")
                    bed_tool = pybedtools.BedTool(bamfile)
                    bedfile_list.append(bed_tool.bam_to_bed.moveto(bamfile.replace("_sorted.bam", ".bed")))

                    # bedfile = bamfile.replace("_sorted.bam", ".bed")
                    # bedfile_list.append(bedfile)
                    # cmd = "bedtools bamtobed -i {0} > {1}" .format(bamfile, bedfile)
                    # os.system(cmd)
                    count += 1
                    print("\t-->Wrote \033[1;214m{0}[m.".format(bamfile.replace("_sorted.bam", ".bed")))

        # I see no reason to keep all the demultiplexed BAM files.
        Tool_Box.delete(self.bamfile_list)

        print("BAM to BED file conversion complete.")

        return bedfile_list
