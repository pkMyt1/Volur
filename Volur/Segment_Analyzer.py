# -*- coding: utf-8 -*-
"""
Segment_Analyzer.py
Version 0.10.0
    February 25, 2019
    Dennis A. Simpson
    Added code to prevent segments that define abrupt CNV changes (del to ins or ins to del) from being counted twice.
Version 0.9.0
    February 17, 2017
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).  Cleaned up code.  Added .format() for
    string building.
Version 0.9
    July 24, 2016
    Dennis A. Simpson
    The permutation analysis is now in its own module.
Version 0.8
    May 8, 2016
    Dennis A. Simpson
    Added ability to do a permutation/frequency type analysis similar to the coin toss module in R.  The intent is to
    use this to calculate the chance occurrence of the observed break point results.  I am not sure that the method
    employed gives the intended results.
Version 0.7
    February 10, 2016
    Dennis A. Simpson
    Aberrations now pick up aneuploid chromosomes.  Changed from Python intersect to numpy intersect for finding common
    aberrations.  Optimized the code.  Reduced runtime on Win10 laptop by about 8 minutes.
Version 0.5
    January 31, 2016
    Dennis A. Simpson
    Program outputs a gene list now.  Format of list and cutoff for percentage of positive cells required needs to be
    modified.  Need to consider binning the gene list by cutoff percentage.
Version 0.2
    January 20, 2016
    Dennis A. Simpson
    Minor tweaks to program structure that altered how the output file name is generated.  Makes this module more
    consistent with how files are named in the other modules.
Version 0.1
    Created on Thu Dec  3 10:03:15 2015
    Dennis A. Simpson
    Takes as input the tab delimited SegCopy file from GINKO and outputs average ploidy per chromosome, size
    distribution and number of insertions and deletions for each chromosome and cell.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2016
"""

import os
import natsort
import statistics
from contextlib import suppress
import Valkyries.Tool_Box as Tool_Box
import Valkyries.Ensembl as Ensembl
from collections import Counter, defaultdict, namedtuple
import numpy
import pybedtools
from datetime import date

__author__ = 'Dennis A. Simpson'
__version__ = "0.10.0"
__package__ = 'Mimir'


class SegmentAnalyzer:

    def __init__(self, log, args):

        self.log = log
        self.args = args
        self.seg_copy_array = None
        self.bin_tracking_array = None
        self.chr_ploidy_dict = defaultdict(lambda: defaultdict(tuple))
        self.abr_range_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.__dict_build()

        if not eval(args.Permutation_Analysis):
            self.ref_data = Ensembl.initialize(args, log)

        # initialize these but only use when Permutation_Analysis calls.
        self.break_point_list = []
        self.unique_permutation_group_size = 0
        self.total_permutation_group_size = 0

    @property
    def chrom_list(self):
        return Tool_Box.chromosomes(self.args.Species, self.log, eval(self.args.Include_chrY))

    @property
    def sample_names(self):
        return self.seg_copy_array[0, 4:].tolist()

    def target_intersection(self, permuted_group=None):
        """
        This will output the intersection of a target bed file and the breakpoints of a cell type.  Output is both by
        breakpoint count and by target count.

        :return:
        """

        def targeting(seg_copy_array, chrom_list, cell_name, abr_range_dict, cell_names, pmuted):
            outstring = ""
            tot_breakpoints = defaultdict(list)
            breakpnt_coordinate_dict = defaultdict(lambda: defaultdict(str))
            seg_counts_dict = defaultdict(lambda: defaultdict(int))

            shuffle_dict_unpaired = defaultdict(list)
            list_length = seg_copy_array.shape[0]
            for cell in sample_names:
                shuffle_dict_unpaired[cell.decode()] = [0] * list_length

            for chrm in natsort.natsorted(chrom_list):
                chrom_slice = seg_copy_array[seg_copy_array[:, 1] == chrm.encode()]

                for cell in natsort.natsorted(cell_names):
                    if not pmuted:
                        with suppress(IndexError):
                            i = len(cell_name)
                            cell_label = cell[:i]
                    else:
                        cell_label = cell
                        cell_name = cell

                    if cell_label == cell_name:
                        with suppress(AttributeError):
                            cell = cell.decode()

                        for aberration in abr_range_dict[chrm][cell]:
                            if abr_range_dict[chrm][cell][aberration]:
                                for segments in abr_range_dict[chrm][cell][aberration]:
                                    for i in range(len(segments)):
                                        process_stop = True

                                        # Prevents us from moving off the 5' end of a chromosome.
                                        try:
                                            coord_start = int(chrom_slice[chrom_slice[:, 0] == (segments[i]) - 1][0, 2])

                                        except IndexError:
                                            coord_start = int(chrom_slice[chrom_slice[:, 0] == segments[i]][0, 2])
                                            coord_stop = int(chrom_slice[chrom_slice[:, 0] == (segments[i]) + 1][0, 3])
                                            process_stop = False

                                        # Prevents us from moving off the 3' end of a chromosome.
                                        if process_stop:
                                            try:
                                                coord_stop = int(chrom_slice[chrom_slice[:, 0] == segments[i]][0, 3])
                                            except IndexError:
                                                coord_stop = int(
                                                    chrom_slice[chrom_slice[:, 0] == (segments[i]) - 1][0, 3])
                                                coord_start = int(
                                                    chrom_slice[chrom_slice[:, 0] == (segments[i]) - 2][0, 2])

                                        try:
                                            # This prevents segments that define a change from del to ins being counted twice.
                                            if not shuffle_dict_unpaired[cell][segments[i]] == 1:
                                                outstring += "{0} {1} {2} {3} {0}|{1}|{2}|{3}\n" \
                                                    .format(chrm, coord_start, coord_stop, aberration)

                                                ky = "{}.{}".format(chrm, coord_start)
                                                seg_counts_dict[ky][aberration] += 1
                                                tot_breakpoints[aberration].append(ky)
                                                breakpnt_coordinate_dict[ky][aberration] += \
                                                    "{0}\t{1}\t{2}\n".format(chrm, coord_start, coord_stop)

                                        except IndexError:
                                            pass

                                        shuffle_dict_unpaired[cell][segments[i]] = 1

            return tot_breakpoints, seg_counts_dict, breakpnt_coordinate_dict, outstring

        if eval(self.args.Target_Intersect_File) and eval(self.args.Breakpoint_Coordinate_File):
            self.log.info("Begin Compiling Breakpoint Coordinates and Breakpoint Target Intersects")
        elif eval(self.args.Target_Intersect_File):
            self.log.info("Begin Compiling Breakpoint Target Intersects")

        if permuted_group is None:
            sample_names = self.sample_names
            encoded_cell_name = self.args.Cell_Name.encode()
            permuted = False
        else:
            sample_names = permuted_group
            encoded_cell_name = None
            permuted = True

        total_breakpoints, seg_count_dict, breakpoint_coordinate_dict, bed_string = \
            targeting(self.seg_copy_array, self.chrom_list, encoded_cell_name, self.abr_range_dict, sample_names,
                      permuted)

        ins_outstring = "#Chrom\tSeg_Start\tSeg_End\n"
        del_outstring = "#Chrom\tSeg_Start\tSeg_End\n"

        if eval(self.args.Breakpoint_Coordinate_File):
            for k in natsort.natsorted(breakpoint_coordinate_dict):
                del_outstring += breakpoint_coordinate_dict[k]["del"]
                ins_outstring += breakpoint_coordinate_dict[k]["ins"]
            cell_name = getattr(self.args, "Cell_Name", "")
            ins_outfile = open("{0}{1}_{2}_ins_breakpoint.bed"
                               .format(self.args.Working_Folder, self.args.Job_Name, cell_name), 'w')
            del_outfile = open("{0}{1}_{2}_del_breakpoint.bed"
                               .format(self.args.Working_Folder, self.args.Job_Name, cell_name), 'w')

            del_outfile.write(del_outstring)
            ins_outfile.write(ins_outstring)
            del_outfile.close()
            ins_outfile.close()

            self.log.info("Breakpoint Coordinate BED Files Written")

        # Only do the pybedtools intersect if the user requests an intersect file or sample permutation.
        if not eval(self.args.Target_Intersect_File) and permuted_group is None:
            return

        self.log.debug("Working on target intersects")

        # Bedtool Section.
        breakpoint_bedtool = pybedtools.BedTool(bed_string, from_string=True)
        target_bedtool = pybedtools.BedTool(self.args.Target_File, from_string=False)

        # Find target intersects for printing.
        breakpoint_target_intersect = breakpoint_bedtool.intersect(target_bedtool, wb=True, stream=True)

        """
        The breakpoint target intersect pybedtools object is expected to have this structure;
        l[0] = Breakpoint chrom; l[1] = Breakpoint start coord; l[2] = Breakpoint end coord; 
        l[3] = aberration copy type; l[4] = segment ID for internal tracking.  The next items are from the target BED 
        file.  Make sure column 5 in that file is the target name.
        """

        intersect_dict = defaultdict(lambda: defaultdict(list))

        # Processing Breakpoint Intersects.
        for l in breakpoint_target_intersect:
            chrom = l[4].split("|")[0]
            start = l[4].split("|")[1]
            stop = l[4].split("|")[2]
            c_type = l[4].split("|")[3]
            k1 = "{}.{}.{}".format(chrom, start, stop)

            intersect_dict[k1][c_type].append(l[9])

        # Format Data for Printing.
        total_targeted_del_breakpoints = 0
        total_targeted_ins_breakpoints = 0
        unique_targeted_del_breakpoints = 0
        unique_targeted_ins_breakpoints = 0
        total_del_targets = 0
        total_ins_targets = 0
        unique_ins_targets = 0
        unique_del_targets = 0

        unique_breakpoint_string = \
            "\nCell\tChrom\tStart\tStop\tCopyType\tTotalSegCount\tTotalTargetCount\tUniqueTargetCount\tTargets\n"

        for key in natsort.natsorted(intersect_dict):
            chrom = key.split(".")[0]
            start = key.split(".")[1]
            stop = key.split(".")[2]
            seg_key = "{}.{}".format(chrom, start)

            total_targeted_del_breakpoints += seg_count_dict[seg_key]["del"]
            if seg_count_dict[seg_key]["del"] > 0:
                unique_targeted_del_breakpoints += 1

            total_targeted_ins_breakpoints += seg_count_dict[seg_key]["ins"]
            if seg_count_dict[seg_key]["ins"] > 0:
                unique_targeted_ins_breakpoints += 1

            total_del_targets += len(intersect_dict[key]["del"])
            total_ins_targets += len(intersect_dict[key]["ins"])
            unique_del_targets += len(set(intersect_dict[key]["del"]))
            unique_ins_targets += len(set(intersect_dict[key]["ins"]))

            for c_type in intersect_dict[key]:
                if len(set(intersect_dict[key][c_type])) > 0:
                    target_string = "\t".join(list(set(intersect_dict[key][c_type])))

                    cell_name = getattr(self.args, "Cell_Name", "")
                    unique_breakpoint_string += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n" \
                        .format(cell_name, chrom, start, stop, c_type, seg_count_dict[seg_key][c_type],
                                len(intersect_dict[key][c_type]),  len(set(intersect_dict[key][c_type])), target_string)

        # Pull the untargeted breakpoints.
        total_del = len(total_breakpoints["del"])
        total_unique_del = len(set(total_breakpoints["del"]))
        total_ins = len(total_breakpoints["ins"])
        total_unique_ins = len(set(total_breakpoints["ins"]))

        if eval(self.args.Permutation_Analysis):
            Breakpoints = \
                namedtuple('Breakpoints', ['total_del', 'total_ins', 'total_targeted_del_breakpoints',
                                           'total_targeted_ins_breakpoints', 'total_unique_del', 'total_unique_ins',
                                           'unique_targeted_del_breakpoints', 'unique_targeted_ins_breakpoints'])

            return Breakpoints(total_del, total_ins, total_targeted_del_breakpoints, total_targeted_ins_breakpoints,
                               total_unique_del, total_unique_ins, unique_targeted_del_breakpoints,
                               unique_targeted_ins_breakpoints)

        breakintersect_outfile = open("{0}{1}_{2}_breakpoint_intersect.txt"
                                      .format(self.args.Working_Folder, self.args.Job_Name, self.args.Cell_Name), 'w')

        breakintersect_string = "Segment Analyzer Target Intersect Module run {0}\n" \
                                "Target File:\t{1}\nSegCopy File:\t{2}\n\n" \
                                "Cell\tCopyType\tTotalBreakpoints\tTotalTargetedBreakpoints\tTotalTargets\t" \
                                "UniqueTotalBreakpoints\tUniqueTargetedBreakpoints\tUniqueTargets\n" \
            .format(date.today().strftime("%Y-%m-%d"), self.args.Target_File, self.args.Segment_File)

        breakintersect_string += \
            "{0}\tdel\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n{0}\tins\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n{13}"\
            .format(self.args.Cell_Name, total_del, total_targeted_del_breakpoints, total_del_targets, total_unique_del,
                    unique_targeted_del_breakpoints, unique_del_targets, total_ins, total_targeted_ins_breakpoints, total_ins_targets,
                    total_unique_ins, unique_targeted_ins_breakpoints, unique_ins_targets, unique_breakpoint_string)

        breakintersect_outfile.write(breakintersect_string)
        breakintersect_outfile.close()
        self.log.info("Breakpoint Target Intersects File Written.")

        return

    def bin_sizing(self):
        """
        This will re-bin the genome space to be consistent with the breakpoint definition or any other bin re-sampling
        desired.
        :return:
        """

        self.log.info("Begin Re-Binning the Genome Space.")
        bed_group = int(self.args.bed_group)
        new_list = []
        seg_num = 0

        for chrom in natsort.natsorted(self.chrom_list):
            self.log.info("Binning Chromosome {0}".format(chrom))

            chrom_slice = self.seg_copy_array[self.seg_copy_array[:, 1] == chrom.encode()]
            seg_count = chrom_slice.shape[0]
            coord_start = int(chrom_slice[0, 2])

            for i in range((seg_count-1)):
                if (i+1) < seg_count and (i+1) % bed_group == 0:
                    coord_stop = int(chrom_slice[i, 3])
                    new_list.append([seg_num, chrom.encode(), coord_start, coord_stop])
                    seg_num += 1

                    # This prevents running past the end of the chromosome.
                    try:
                        coord_start = int(chrom_slice[i+bed_group, 2])
                    except IndexError:
                        coord_start = int(chrom_slice[seg_count-1, 2])

        self.log.info("Genome Space Successfully Re-Binned.")

        return numpy.array(new_list, dtype='object')

    @staticmethod
    def sub_target_mapping(permutation_data, chrom):
        """
        This is the subprocess that computes the mapping coordinates for a target bed file on a single chromosome.
        :param permutation_data:
        :param chrom:
        :return:
        """

        print("\tGenerating Map Key for Chromosome {0}".format(chrom))
        Tool_Box.debug_messenger("Is this ever used?")

        map_list = []
        bin_tracking_array = permutation_data[0]
        target_bed_array = permutation_data[1]

        chrom_slice = bin_tracking_array[bin_tracking_array[:, 1] == chrom.encode()]
        seg_count = chrom_slice.shape[0]
        target_slice = target_bed_array[target_bed_array[:, 1] == chrom.encode()]
        target_row_count = target_slice.shape[0]

        for i in range(seg_count-1):
            coord_start = int(chrom_slice[i, 2])
            coord_stop = int(chrom_slice[i, 3])
            source_seg_id = int(chrom_slice[i, 0])

            target_id_list = []
            seg_match = False

            for j in range(target_row_count-1):
                target_start = int(target_slice[j, 2])
                target_stop = int(target_slice[j, 3])
                target_seg_id = int(target_slice[j, 0])

                if coord_start <= target_start <= coord_stop:
                    seg_match = True
                    target_id_list.append(target_seg_id)

                elif coord_stop >= target_stop >= coord_start:
                    seg_match = True
                    target_id_list.append(target_seg_id)

                elif coord_start >= target_stop >= coord_stop:
                    seg_match = True
                    target_id_list.append(target_seg_id)

            if seg_match:
                map_list.append([source_seg_id, tuple(target_id_list)])

        return map_list

    def __dict_build(self):
        """
        Building several dictionaries here.
        :return:
        """

        self.log.debug("Initializing Dictionaries and Arrays")

        if not os.path.isfile(self.args.Segment_File):
            self.log.error("Segment_File Not Found.  Check File Name and Path.")
            raise SystemExit(1)

        # Parse segment file into array.
        self.log.debug("Parse Segment File into Array")
        initial_array = numpy.genfromtxt(self.args.Segment_File, delimiter='\t', dtype=object)

        array_index_list = []
        for i in range(initial_array.shape[0]):
            array_index_list.append([i])

        # Append index column to initial array.
        self.seg_copy_array = numpy.append(array_index_list, numpy.asarray(initial_array), axis=1)

        array_index_list.clear()
        tmp_list = []

        for row in self.seg_copy_array:
            tmp_list.append([row[0], row[1], row[2], row[3]])

        self.bin_tracking_array = numpy.array(tmp_list, dtype='object')

        self.log.debug("Segment File Array Created")
        self.log.debug("Dictionaries and Arrays Initialized")

    def chromosome_ploidy(self, args=None, permutation=False):
        """
        This will calculate the ploidy for each chromosome and write it to an output file.  Each row is a cell,
        the columns are the ploidy for the cell and then the ploidy for each chromosome analyzed.  Functions that
        require the chr_ploidy_dict are dependent on this function.

        :return:
        """

        self.log.info("Begin Calculating Chromosome Ploidy.")

        if args:
            self.args = args

        tmp_ploidy_dict = defaultdict(lambda: defaultdict(list))
        ploidy_out_dict = defaultdict(list)
        cell_labels = self.seg_copy_array[self.seg_copy_array[:, 0] == 0]

        # Read data from segment copy array one chromosome at a time.
        for chrom in natsort.natsorted(self.chrom_list):
            self.log.debug("Analyzing {0}".format(chrom))
            chrom_slice = self.seg_copy_array[self.seg_copy_array[:, 1] == chrom.encode()]
            row_count = self.seg_copy_array[self.seg_copy_array[:, 1] == chrom.encode()].shape[0]
            column_count = self.seg_copy_array[self.seg_copy_array[:, 1] == chrom.encode()].shape[1]

            for r in range(row_count):
                for i in range(4, column_count):
                    cell = cell_labels[0, i].decode()
                    ploidy = int(chrom_slice[r, i])
                    tmp_ploidy_dict[chrom][cell].append(ploidy)

        # Use data in temporary dictionary to calculate the ploidy of each chromosome for each cell.
        self.log.info("Begin Processing Ploidy Stats.")
        for chrom in natsort.natsorted(tmp_ploidy_dict):
            self.log.debug("Processing Ploidy for {}".format(chrom))
            for cell, v in natsort.natsorted(tmp_ploidy_dict[chrom].items()):
                avg = statistics.mean(v)
                mdian = statistics.median(v)
                std = statistics.pstdev(v)
                self.chr_ploidy_dict[chrom][cell] = avg
                ploidy_out_dict[cell].append((chrom, avg, std, mdian))

        if not permutation and eval(self.args.Chrom_Ploidy_File):
            self.log.info("Writing {}_ploidy_distribution.txt".format(self.args.Job_Name))
            if eval(self.args.Include_chrY):
                h2 = "\tChrY_Mean\tStDev\tMedian\n"
            else:
                h2 = "\n"

            ploidy_outstring = ("Cell\tCell_Ploidy"
                                "\tChr1_Mean\tStDev\tMedian\tChr2_Mean\tStDev\tMedian\tChr3_Mean\tStDev\tMedian"
                                "\tChr4_Mean\tStDev\tMedian\tChr5_Mean\tStDev\tMedian\tChr6_Mean\tStDev\tMedian"
                                "\tChr7_Mean\tStDev\tMedian\tChr8_Mean\tStDev\tMedian\tChr9_Mean\tStDev\tMedian"
                                "\tChr10_Mean\tStDev\tMedian\tChr11_Mean\tStDev\tMedian\tChr12_Mean\tStDev\tMedian"
                                "\tChr13_Mean\tStDev\tMedian\tChr14_Mean\tStDev\tMedian\tChr15_Mean\tStDev\tMedian"
                                "\tChr16_Mean\tStDev\tMedian\tChr17_Mean\tStDev\tMedian\tChr18_Mean\tStDev\tMedian"
                                "\tChr19_Mean\tStDev\tMedian\tChrX_Mean\tStDev\tMedian{0}".format(h2))

            # Iterate the ploidy dictionary and write the data output file.
            for cell in natsort.natsorted(ploidy_out_dict):
                param_list = ""
                ploidy_sum = 0
                chrom_count = 0
                for ploidy_list in natsort.natsorted(ploidy_out_dict[cell]):
                    chrom_count += 1
                    ploidy_sum += ploidy_list[1]
                    param_list += "{0}\t{1}\t{2}\t".format(round(ploidy_list[1], 2), round(ploidy_list[2], 2),
                                                           ploidy_list[3])
                cell_ploidy = int(round(ploidy_sum/chrom_count, 0))
                ploidy_outstring += "{0}\t{1}\t{2}\n".format(cell, cell_ploidy, param_list)

            ploidy_outfile = \
                open("{0}{1}_ploidy_distribution.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
            ploidy_outfile.write(ploidy_outstring)
            ploidy_outfile.close()
            self.log.info("{}_ploidy_distribution.txt Written".format(self.args.Job_Name))

        # Clear dictionaries to free up memory.
        tmp_ploidy_dict.clear()
        ploidy_out_dict.clear()

        self.log.info("Chromosome Ploidy Analysis Complete.")

        return

    def break_points(self, cell=None, permutation=False):
        """
        This function will process the breakpoints and return data files by chromosome or segment.  The segments contain
        the bin prior to the breakpoint as well.  Also used by Permutation_Analysis to get the number of breakpoints per
        cell type.
        :param cell:
        :param permutation:
        :return:
        """

        if not permutation:
            self.log.info("Begin Break Point Calculations.")
        else:
            label = ""
            if self.args.Permutation_Type == "Segment":
                self.args.Cell_Name = cell
                label = " on {}".format(self.args.Cell_Name)

            self.log.info("Segment Analyzer Determining Breakpoints for {} Type Permutation Analysis{}"
                          .format(self.args.Permutation_Type, label))

        # Initialize some of our variables.
        param_string = ""
        bp_chrom_distribution_data = ""
        seg_copy_array = self.seg_copy_array
        cell_labels = seg_copy_array[seg_copy_array[:, 0] == 0]

        # Begin the analysis.  Proceeds one chromosome at a time.
        for chrom in natsort.natsorted(self.chrom_list):
            if not permutation:
                self.log.debug("Begin Compiling Break Points on {0} for {1}".format(chrom, self.args.Cell_Name))

            # Copy Type and Segment Size calls.
            tmp_seg_dict, cell_break_dict = self.copy_type(chrom, seg_copy_array, cell_labels)

            # Iterate the segment dictionary.  Make a breakpoint list, count number of unique breakpoints, format data
            # for segment breakpoint file.
            for seg in sorted(tmp_seg_dict):
                # ToDo: Clean this block up.
                self.break_point_list.append(seg)
                self.unique_permutation_group_size += 1
                for abr_type in tmp_seg_dict[seg]:
                    self.total_permutation_group_size += tmp_seg_dict[seg][abr_type]

                # This block formats the data for the segmented breakpoint file.
                if not permutation:
                    if eval(self.args.Breakpoint_Dist_File):
                        for abr_type in sorted(tmp_seg_dict[seg]):
                            counts = tmp_seg_dict[seg][abr_type]

                            if counts > 0:
                                coord_start = int(self.seg_copy_array[self.seg_copy_array[:, 0] == (seg-1)][0, 2])
                                coord_stop = int(self.seg_copy_array[self.seg_copy_array[:, 0] == seg][0, 3])

                                param_string += "{0}\t{1}\t{2}\t{3}\t{4}\n"\
                                    .format(chrom, coord_start, coord_stop, counts, abr_type)

            # This block formats the data for the chromosome breakpoints file.
            if not permutation:
                if eval(self.args.Breakpoint_Chrom_Dist_File):
                    for cell_name in natsort.natsorted(cell_break_dict):
                        del_count = cell_break_dict[cell_name]["del"]
                        ins_count = cell_break_dict[cell_name]["ins"]

                        bp_chrom_distribution_data += "{0}\t{1}\t{2}\t{3}\t{4}\n"\
                            .format(chrom, cell_name, del_count, ins_count, del_count+ins_count)

        # Write the data to the appropriate output files.
        if not permutation:
            if eval(self.args.Breakpoint_Dist_File):
                file_outstring = "#chrom\tstart\tstop\tcount\tabr_type\n"
                breakpoint_distribution = open("{0}{1}_{2}_segmented_breakPoint.bed"
                                               .format(self.args.Working_Folder, self.args.Job_Name,
                                                       self.args.Cell_Name), 'w')

                file_outstring += param_string
                breakpoint_distribution.write(file_outstring)
                breakpoint_distribution.close()

            if eval(self.args.Breakpoint_Chrom_Dist_File):
                file_append2 = "allCell_chromosome_breakPoint.bed"
                file_outstring = "#chrom\tcell_name\tdel_count\tins_count\ttotal\n"
                bp_chrom_distribution = \
                    open("{0}{1}_{2}".format(self.args.Working_Folder, self.args.Job_Name, file_append2), 'w')
                file_outstring += bp_chrom_distribution_data
                bp_chrom_distribution.write(file_outstring)
                bp_chrom_distribution.close()

        if not permutation:
            self.log.info("Break Point Calculations Complete.")

        return

    def copy_type(self, chrom, seg_copy_array, cell_labels):

        chrom_slice = seg_copy_array[seg_copy_array[:, 1] == chrom.encode()]
        row_count = chrom_slice.shape[0]
        column_count = chrom_slice.shape[1]
        prior_ctype_dict = defaultdict(str)
        prior_ploidy = defaultdict(list)
        tmp_seg_dict = defaultdict(lambda: defaultdict(int))
        cell_break_dict = defaultdict(lambda: defaultdict(int))

        for r in range(row_count):
            bin_num = int(chrom_slice[r, 0])
            # bin_start = int(chrom_slice[r, 2])

            for c in range(4, column_count):

                cell = cell_labels[0, c].decode()
                cell_label = cell
                cell_name = cell

                if eval(self.args.Permutation_Analysis):
                    if self.args.Permutation_Type == "Segment":
                        with suppress(IndexError):
                            i = len(self.args.Cell_Name)
                            cell_label = cell[:i]
                            cell_name = self.args.Cell_Name

                    elif not self.args.Permutation_Type == "Sample" \
                            and not self.args.Permutation_Type == "Segment" \
                            and not self.args.Permutation_Type == "Shuffle":

                        self.log.error("--Permutation_Type option not Sample, Segment or Shuffle.")
                        raise SystemExit(1)

                # Segment Permutations are done on single cell types.
                if not cell_label == cell_name:
                    continue

                current_ploidy = int(chrom_slice[r, c])  # Set current ploidy equal to segment ploidy.

                # Define the current segment as normal, insertion or deletion.
                if current_ploidy == round(self.chr_ploidy_dict[chrom][cell]):
                    copy_type = ""
                elif current_ploidy < round(self.chr_ploidy_dict[chrom][cell]):
                    copy_type = "del"
                elif current_ploidy > round(self.chr_ploidy_dict[chrom][cell]):
                    copy_type = "ins"

                # If the copy type changes or we reach the end of the chromosome with abnormal copy type we
                # need to do some things.
                if not copy_type == prior_ctype_dict[cell] or not copy_type == "" and r == row_count - 1:
                    # Flags are used when copy type goes from insertion to deletion without a chrom average.
                    flag_ins = False
                    flag_del = False

                    if r == 0:
                        # First row for a chromosome.
                        c_type = copy_type
                    elif copy_type == "":
                        c_type = prior_ctype_dict[cell]
                    elif copy_type == "ins" and prior_ctype_dict[cell] == "":
                        # New insertion
                        c_type = copy_type
                    elif copy_type == "ins" and prior_ctype_dict[cell] == "del":
                        # Deletion going to insertion
                        flag_ins = True
                        c_type = prior_ctype_dict[cell]
                    elif copy_type == "del" and prior_ctype_dict[cell] == "":
                        # New Deletion
                        c_type = copy_type
                    elif copy_type == "del" and prior_ctype_dict[cell] == "ins":
                        # Insertion going to deletion
                        flag_del = True
                        c_type = prior_ctype_dict[cell]
                    elif r == row_count - 1:
                        # End of chromosome
                        c_type = copy_type
                        prior_ploidy[cell][1] = True
                    else:
                        self.log.eror(chrom, cell, chrom_slice[r, 0], copy_type, current_ploidy,
                                      round(self.chr_ploidy_dict[chrom][cell]))
                        Tool_Box.debug_messenger("There is an error in the segment copy file.")

                        raise SystemExit(1)

                    # Build the segment breakpoint dictionary and load data into both the segment and chromosome
                    # dictionaries.
                    seg = chrom_slice[r, 0]

                    if eval(self.args.Permutation_Analysis) and self.args.Permutation_Type == "Sample":
                        # Process data for Sample type Permutation Analysis.
                        if seg == chrom_slice[0, 0]:
                            s1 = [seg, seg + 1]
                        else:
                            s1 = [seg - 1, seg]

                        self.break_point_list.extend(s1)

                        if flag_del or flag_ins:
                            s1 = [seg + 1, seg + 2]
                            self.break_point_list.extend(s1)
                    else:
                        # Process data for Segment Analysis or Segment type Permutation Analysis.
                        if r == 0:
                            # There is no previous segment if at the beginning so step forward 1.
                            seg += 1

                        if cell[:len(self.args.Cell_Name)] == self.args.Cell_Name:
                            tmp_seg_dict[seg][c_type] += 1

                        cell_break_dict[cell_name][c_type] += 1

                        if flag_del or flag_ins:
                            if r == row_count - 1:
                                seg -= 1

                            if flag_ins:
                                c_type = "ins"
                            elif flag_del:
                                c_type = "del"

                            cell_break_dict[cell_name][c_type] += 1
                            if cell[:len(self.args.Cell_Name)] == self.args.Cell_Name:
                                tmp_seg_dict[seg + 1][c_type] += 1

                prior_ctype_dict[cell] = copy_type

                # Do some Segment Size Measurements
                if r == 0:
                    prior_ploidy[cell] = [copy_type, False, bin_num]

                elif not copy_type == prior_ploidy[cell][0]:
                    if prior_ploidy[cell][0] == "":
                        prior_ploidy[cell] = [copy_type, False, bin_num]
                    else:
                        prior_ploidy[cell][1] = True

                elif copy_type == "ins" and prior_ploidy[cell][0] == "del":
                    prior_ploidy[cell][0] = "del"
                    prior_ploidy[cell][1] = True

                elif copy_type == "del" and prior_ploidy[cell][0] == "ins":
                    prior_ploidy[cell][0] = "ins"
                    prior_ploidy[cell][1] = True

                elif r == row_count - 1 and (copy_type == 'del' or copy_type == "ins"):
                    # End of chromosome
                    prior_ploidy[cell][1] = True

                if prior_ploidy[cell][1]:
                    # Add the data to the dictionaries.
                    self.abr_range_dict[chrom][cell][prior_ploidy[cell][0]].extend([(prior_ploidy[cell][2], bin_num)])
                    prior_ploidy[cell] = [copy_type, False, bin_num]

        return tmp_seg_dict, cell_break_dict

    def genes(self):
        """
        This generates a list of genes overlapping the CNV region.
        :param self:
        :return:
        """
        self.log.info("Begin Gene Mapping")

        outstring = "Recurrent Count:\t≥{0}\nCell Type:\t{1}\nChrom\tAberration\tGenes\n"\
            .format(self.args.Minimum_Family_Size, self.args.Cell_Name)

        # ToDo: Put this into Pandas so the columns will be chromosomes and the rows the genes.
        for chrom in natsort.natsorted(self.chrom_list):
            chrom_slice = self.seg_copy_array[self.seg_copy_array[:, 1] == chrom.encode()]
            full_gene_aberration_dict = defaultdict(list)

            for cell in natsort.natsorted(self.abr_range_dict[chrom]):

                with suppress(IndexError):
                    i = len(self.args.Cell_Name)
                    cell_label = cell[:i]

                # Restrict to cell type from options file.  Generate a list of genes for each of these cells.
                if cell_label == self.args.Cell_Name:
                    for aberration in self.abr_range_dict[chrom][cell]:
                        if self.abr_range_dict[chrom][cell][aberration]:
                            for segments in self.abr_range_dict[chrom][cell][aberration]:

                                start = int(chrom_slice[chrom_slice[:, 0] == segments[0]][0, 2])
                                stop = int(chrom_slice[chrom_slice[:, 0] == segments[1]][0, 3])
                                full_gene_aberration_dict[aberration].extend(
                                    self.ref_data.gene_names_at_locus(contig=chrom, position=start, end=stop))

            # Find the recurrent genes.
            outstring = \
                self._recurrent_target_formatting(chrom, Counter(full_gene_aberration_dict["del"]), "del", outstring)
            outstring = \
                self._recurrent_target_formatting(chrom, Counter(full_gene_aberration_dict["ins"]), "ins", outstring)

        recurrent_target_outfile = \
            open("{0}{1}_{2}_RCount{3}_targets.txt"
                 .format(self.args.Working_Folder, self.args.Job_Name, self.args.Cell_Name,
                         self.args.Minimum_Family_Size), 'w')

        recurrent_target_outfile.write(outstring)
        recurrent_target_outfile.close()

        return

    def _recurrent_target_formatting(self, chrom, counter_obj, c_type, outstring):
        results = False
        gene_str = ""
        for gene in counter_obj:
            if counter_obj[gene] >= int(self.args.Minimum_Family_Size):
                results = True
                gene_str += "{0}\t".format(gene)

        if results:
            outstring += "{0}\t{1}\t{2}\n".format(chrom, c_type, gene_str)

        return outstring

    def size_output(self):
        """
        This will generate a file of CNV's with location and type.
        :return:
        """

        self.log.info("Begin Size Distribution Analysis")
        size_distribution = open("{0}{1}_size_distribution.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        size_distribution.write("Chrom\tCell\tType\tCount\tTotal_Size(KB)\tSize(KB)")
        outstring = ""

        # Calculate the size of the aberrations in KB and build a dictionary of the results.
        for chrom in natsort.natsorted(self.abr_range_dict):
            self.log.debug("Processing {0} for aberration size output.".format(chrom))

            for cell in natsort.natsorted(self.abr_range_dict[chrom]):
                tmp_outstring = ""
                print_flag = False

                for abr_type in natsort.natsorted(self.abr_range_dict[chrom][cell]):
                    abr_count = len(self.abr_range_dict[chrom][cell][abr_type])
                    abr_list = []

                    if abr_count > 0:
                        print_flag = True

                    for segments in self.abr_range_dict[chrom][cell][abr_type]:
                        coord_start = int(self.seg_copy_array[self.seg_copy_array[:, 0] == segments[0]][0, 2])
                        chrom_start = self.seg_copy_array[self.seg_copy_array[:, 0] == segments[0]][0, 1].decode()
                        coord_stop = int(self.seg_copy_array[self.seg_copy_array[:, 0] == segments[1]][0, 3])
                        chrom_stop = self.seg_copy_array[self.seg_copy_array[:, 0] == segments[1]][0, 1].decode()
                        if not chrom_start == chrom or not chrom_stop == chrom:
                            self.log.error("Segment File error at segments {}, chr{}:{}-{}"
                                           .format(segments, chrom, coord_start, coord_stop))
                            raise SystemExit("1")

                        abr_list.append(int((coord_stop - coord_start) / 1000))

                    abr_list_str = "\t".join(str(value) for value in sorted(abr_list))
                    abr_size = sum(abr_list)
                    tmp_outstring += "\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}"\
                        .format(chrom, cell, abr_type, abr_count, abr_size, abr_list_str)

                if print_flag:
                    outstring += tmp_outstring

        size_distribution.write(outstring)
        size_distribution.close()

        self.log.info("Size Distribution Analysis Complete.  {}_size_distribution.txt written"
                      .format(self.args.Job_Name))

        return
