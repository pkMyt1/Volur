"""
Version 0.10.0
    March 1, 2019
    Dennis A. Simpson
    Added the ability to keep the breakpoint pairs together when doing shuffle permutations.  Requires
    --PairedBreakpoints/tTrue or False in the parameter file.
Version 0.9.0
    February 22, 2019
    Dennis A. Simpson
    Added Shuffle Analysis method.  This will shuffle the location of the breakpoints for each cell to create the null
    set.  When it is done it should also be able to keep the paired breakpoint together during the shuffle permutations.
Version 0.8.0
    September 4, 2018
    Dennis A. Simpson
    Split out from Volundr.  Now called Volur.  Only single cell analysis data.  Changed the segment permutation output 
    to include the odds ratio and p-values from Fisher Exact test.
Version 0.5.o
    December 2, 2016
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).  Added version dependency checks for
    modules.
Version 0.5
    July 7, 2016
    Dennis A. Simpson
    Added Target Permutations where the sample space is the target space and the intersection is done with the
    breakpoint bed file.  Also added the ability to control the outputs.  Can now choose any combination of mapfile,
    segment permutation, and target permutation.

Version 0.1
    Created on June 8, 2016
    Dennis A. Simpson
    This class will do the permutation analysis for single cell sequencing.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""

import itertools
from collections import Counter, defaultdict
from contextlib import suppress
import numpy
import natsort
import pathos
from operator import itemgetter
import Volur.Segment_Analyzer as Segment_Analyzer
import Valkyries.Tool_Box as Tool_Box
import scipy
import scipy.stats
import pybedtools
from datetime import date

__author__ = 'Dennis A. Simpson'
__version__ = "0.10.0"
__package__ = 'Volur'


class PermutationAnalysis:
    def __init__(self, log, args):
        if not getattr(args, "Segment_Permutation_File", False) and not getattr(args, "Map_File", False):
            if args.Permutation_Type == "Segment":
                log.error("You have not selected any output file(s) for the Segment Permutation Analysis.")
                raise SystemExit(1)

        self.log = log
        self.args = args
        self.seg_analyzer = self.segment_analyzer()
        self.chromosome_list = Tool_Box.chromosomes(self.args.Species, self.log, eval(self.args.Include_chrY))
        self.shuffle_dict_unpaired = defaultdict(list)
        self.shuffle_dict_pairs = defaultdict(list)
        # initialize but only used if doing Segment Permutation Analysis.
        self.break_points = ""
        self.target_bed_array = None
        self.bin_tracking_array = None
        self.array_build()
        self.target_bed_map_array = None
        self.break_point_list = []
        self.cell_breakpoint_dict = {}
        self.cell_label_list = []

    def shuffle_analysis(self):
        self.log.debug("Calculate breakpoints, and fill aberration dictionary from Segment Analyzer")
        self.seg_analyzer.break_points(cell=None, permutation=True)

        self.log.debug("Populate Cell Label List")
        for name in self.seg_analyzer.sample_names:
            self.cell_label_list.append(name.decode())

        self.log.debug("Begin building shuffle dictionaries")
        list_length = self.seg_analyzer.seg_copy_array.shape[0]
        for cell in self.cell_label_list:
            self.shuffle_dict_pairs[cell] = [0]*list_length
            self.shuffle_dict_unpaired[cell] = [0]*list_length

        self.log.debug("Begin sorting breakpoints into shuffle dictionaries")
        for chrom in natsort.natsorted(self.seg_analyzer.chrom_list):
            self.log.debug("Processing {}".format(chrom))
            # chrom_slice = self.seg_analyzer.seg_copy_array[self.seg_analyzer.seg_copy_array[:, 1] == chrom.encode()]

            for cell in self.cell_label_list:
                for aberration in self.seg_analyzer.abr_range_dict[chrom][cell]:
                    for segments in self.seg_analyzer.abr_range_dict[chrom][cell][aberration]:

                        # Paired Breakpoints
                        self.shuffle_dict_pairs[cell][segments[0]] = segments[1] - segments[0]

                        for i in range(len(segments)):
                            # Unpaired breakpoints
                            self.shuffle_dict_unpaired[cell][segments[i]] = 1

        self.log.info("Spawning {0} jobs to process {1} iterations shuffle permutation analysis."
                      .format(self.args.Spawn, self.args.Iteration_Count))

        idlist = []
        for i in range(int(self.args.Spawn)):
            idlist.append(i)

        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        p.starmap(self.target_intersection, zip(itertools.repeat(self), idlist))

        self.log.info("Shuffle Permutation jobs done.  Compile any temporary data files into final results.")

        outstring = "Shuffle Permuted Analysis for {} Cell Type\n" \
                    "Total Breakpoints\tTotal Intersected\tUnique Breakpoints\tUnique Intersected\n" \
            .format(self.args.Cell_Name)
        file_delete_list = []
        for i in idlist:
            temp_file = "{}{}{}{}".format(self.args.Working_Folder, self.args.Cell_Name, self.args.Job_Name, i)
            file_delete_list.append(temp_file)
            with open(temp_file) as f:
                for l in f:
                    outstring += l
        if eval(self.args.PairedBreakpoints):
            permuted_shuffle_file_name = "{0}{1}_{2}_Paired_Shuffle_Permutation.bed"\
                .format(self.args.Working_Folder, self.args.Job_Name, self.args.Cell_Name)
        else:
            permuted_shuffle_file_name = "{0}{1}_{2}_Unpaired_Shuffle_Permutation.bed"\
                .format(self.args.Working_Folder, self.args.Job_Name, self.args.Cell_Name)
        permuted_shuffle_file = open(permuted_shuffle_file_name, 'w')
        permuted_shuffle_file.write(outstring)
        permuted_shuffle_file.close()
        Tool_Box.delete(file_delete_list)

    @ staticmethod
    def target_intersection(self, runid):
        """
        This will output determine the intersection of a target bed file and the breakpoints of a cell type.  Done as
        part of the Shuffle Permutation Analysis.
        :return:
        """

        def targeting(shuffledict, seg_copy_array, cell_name):
            bedstring = ""
            seg_counts_dict = defaultdict(int)
            breakpoint_counts = 0
            sum_counts = 0
            for cell in shuffledict:
                with suppress(IndexError):
                    i = len(cell_name)
                    cell_label = cell[:i]

                if not cell_name == cell_label:
                    continue

                shuffled_list = shuffledict[cell]
                scipy.random.shuffle(shuffled_list)
                sum_counts += sum(shuffled_list)

                for i in range(len(shuffled_list)):
                    if shuffled_list[i] == 0:
                        continue

                    breakpoint_counts += 1
                    segment_index = i
                    if i == 0:
                        segment_index = 1

                    chrm = seg_copy_array[seg_copy_array[:, 0] == segment_index][0, 1].decode()
                    chrom_slice = seg_copy_array[seg_copy_array[:, 1] == chrm.encode()]
                    chrom_seg_count = chrom_slice.shape[0]
                    start_seg = segment_index
                    stop_seg = segment_index+1

                    # Prevents us from running past the end of the chromosome
                    if segment_index+1 > chrom_seg_count:
                        stop_seg = segment_index
                        start_seg = segment_index-1

                    coord_start = int(seg_copy_array[seg_copy_array[:, 0] == start_seg][0, 2])
                    coord_stop = int(seg_copy_array[seg_copy_array[:, 0] == stop_seg][0, 3])

                    segkey = "{}.{}".format(chrm, coord_start)
                    seg_counts_dict[segkey] += 1
                    bedstring += "{0} {1} {2} {3} {0}|{1}|{2}|{3}\n".format(chrm, coord_start, coord_stop, "x")

                    if eval(self.args.PairedBreakpoints):
                        segment_index = shuffled_list[i]+i

                        # Since segments are paired we can run past the end of the list.
                        if segment_index > len(shuffled_list):
                            segment_index = len(shuffled_list)-1

                        # If the shuffle results in a segment overlap, skip it.
                        if not shuffled_list[segment_index] == 0:
                            continue

                        start_seg = segment_index
                        stop_seg = segment_index+1

                        # Prevents us from running past the end of the chromosome by flipping direction of region
                        if segment_index + 1 > chrom_seg_count:
                            start_seg = shuffled_list[i] - i
                            stop_seg = start_seg-1

                        coor_start = int(seg_copy_array[seg_copy_array[:, 0] == start_seg][0, 2])
                        coor_stop = int(seg_copy_array[seg_copy_array[:, 0] == stop_seg][0, 3])
                        breakpoint_counts += 1
                        segkey = "{}.{}".format(chrm, coord_start)
                        seg_counts_dict[segkey] += 1
                        bedstring += "{0} {1} {2} {3} {0}|{1}|{2}|{3}\n".format(chrm, coor_start, coor_stop, "x")

            return bedstring, seg_counts_dict, breakpoint_counts

        encoded_cell_name = self.args.Cell_Name
        shuffle_dict = self.shuffle_dict_unpaired
        if eval(self.args.PairedBreakpoints):
            shuffle_dict = self.shuffle_dict_pairs
        output_data_dict = defaultdict(lambda: defaultdict(str))

        iteration_limit = int(self.args.Iteration_Count)/int(self.args.Spawn)
        iteration_count = 0
        while iteration_count < iteration_limit:
            if iteration_count % int(self.args.Prog_Check) == 0:
                self.log.info("Iteration: {} of {} for job {}".format(iteration_count, iteration_limit, runid))

            bed_string, segment_count_dict, total_breakpoints = \
                targeting(shuffle_dict, self.seg_analyzer.seg_copy_array, encoded_cell_name)

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

            # Processing Breakpoint Intersects.
            intersect_dict = defaultdict(list)
            total_targeted_breakpoints = 0
            unique_targeted_breakpoints = 0

            for l in breakpoint_target_intersect:
                chrom = l[4].split("|")[0]
                start = l[4].split("|")[1]
                segment_key = "{}.{}".format(chrom, start)
                intersect_dict[segment_key].append(l[9])

            for k in intersect_dict:
                total_targeted_breakpoints += segment_count_dict[k]
                if segment_count_dict[k] > 0:
                    unique_targeted_breakpoints += 1

            output_data_dict[iteration_count] = "{}\t{}\t{}\t{}\n"\
                .format(total_breakpoints, total_targeted_breakpoints, len(segment_count_dict), len(intersect_dict))

            iteration_count += 1

        # Process data for output and write file.
        outstring = ""

        for k in output_data_dict:
            outstring += output_data_dict[k]

        permuted_shuffle_file_name = \
            "{}{}{}{}".format(self.args.Working_Folder, self.args.Cell_Name, self.args.Job_Name, runid)
        permuted_shuffle_file = open(permuted_shuffle_file_name, 'w')
        permuted_shuffle_file.write(outstring)
        permuted_shuffle_file.close()

        return

    def array_build(self):
        self.log.info("Parse Target BED File into Array")
        initial_array = numpy.genfromtxt(self.args.Target_File, delimiter='\t', dtype=object)

        array_index_list = []
        for i in range(initial_array.shape[0]):
            array_index_list.append([i])

        # Append index column to initial array.
        self.target_bed_array = numpy.append(array_index_list, numpy.asarray(initial_array), axis=1)

        array_index_list.clear()

        if int(self.args.Combine_Segments) > 1:
            self.bin_tracking_array = self.bin_sizing()
        elif int(self.args.Combine_Segments) < 1:
            self.log.error("The --Combine_Segments parameter must be >= 1.")
            raise SystemExit(1)

    def segment_analyzer(self):
        """
        Initialize our Segment_Analyzer class to get segment copy array, targeting array, breakpoint counts ect.
        :return:
        """
        sa = Segment_Analyzer.SegmentAnalyzer(self.log, self.args)
        sa.chromosome_ploidy(permutation=True)

        return sa

    def cell_permutation(self):
        """
        Do a permutation analysis based on the cells.  Uses chromosome ploidy and breakpoint modules in Segment_Analyzer
        :return:
        """

        self.log.info("Begin Sample Permutation Analysis.")

        # Initialize some variables.
        self.seg_analyzer.break_points(permutation=True)
        permutation_list = self.seg_analyzer.sample_names
        # cell_permutation_data_dict = defaultdict(lambda: defaultdict(list))
        odds_string = ""
        unique_targeted_odds_ratio_list = []
        total_targeted_odds_ratio_list = []
        total_targeted_del_odds_ratio_list = []
        total_targeted_ins_odds_ratio_list = []
        unique_targeted_ins_odds_ratio_list = []
        unique_targeted_del_odds_ratio_list = []

        # Run a loop for the iterations. Shuffle the list and make a copy for each loop.

        for i in range(int(self.args.Iteration_Count)):
            numpy.random.shuffle(permutation_list)
            shuffled_permutation_list = permutation_list
            sub_list = []
            count = 0

            if i % int(self.args.Prog_Check) == 0:
                self.log.info("Iteration {0} of {1} for Sample Permutation Analysis."
                              .format(i, self.args.Iteration_Count))

                # Pybedtools keeps all temporary files until Python exits.  This helps keep the disk clean.
                pybedtools.cleanup()

            # Create a list with two unique, random lists of indices.
            while count < 2:
                n = (numpy.random.choice(shuffled_permutation_list, int(self.args.Sample_Group_Size), replace=False))

                # Remove the first set from the list
                shuffled_permutation_list = list(set(shuffled_permutation_list).difference(n))
                sub_list.append(n)
                count += 1

            # Retrieve a namedtuple of the permuted samples
            d0 = self.seg_analyzer.target_intersection(sub_list[0])
            d1 = self.seg_analyzer.target_intersection(sub_list[1])

            # cell_permutation_data_dict[0]['del'].append([d0.total_del, d0.total_targeted_del_breakpoints,
            #                                              d0.total_unique_del, d0.unique_targeted_del_breakpoints])
            # cell_permutation_data_dict[1]['del'].append([d1.total_del, d1.total_targeted_del_breakpoints,
            #                                              d1.total_unique_del, d1.unique_targeted_del_breakpoints])
            # cell_permutation_data_dict[0]['ins'].append([d0.total_ins, d0.total_targeted_ins_breakpoints,
            #                                              d0.total_unique_ins, d0.unique_targeted_ins_breakpoints])
            #
            # cell_permutation_data_dict[1]['ins'].append([d1.total_ins, d1.total_targeted_ins_breakpoints,
            #                                              d1.total_unique_ins, d1.unique_targeted_ins_breakpoints])

            total_breakpoint0 = d0.total_del+d0.total_ins
            total_targeted0 = d0.total_targeted_del_breakpoints+d0.total_targeted_ins_breakpoints
            total_unique_breakpoint0 = d0.total_unique_del+d0.total_unique_ins
            total_unique_targeted0 = d0.unique_targeted_del_breakpoints+d0.unique_targeted_ins_breakpoints

            total_breakpoint1 = d1.total_del+d1.total_ins
            total_targeted1 = d1.total_targeted_del_breakpoints+d1.total_targeted_ins_breakpoints
            total_unique_breakpoint1 = d1.total_unique_del+d1.total_unique_ins
            total_unique_targeted1 = d1.unique_targeted_del_breakpoints+d1.unique_targeted_ins_breakpoints

            total_target_ratio0 = total_targeted0/total_breakpoint0
            total_target_ratio1 = total_targeted1/total_breakpoint1

            total_target_odds = total_target_ratio0/total_target_ratio1

            unique_target0 = total_unique_targeted0/total_unique_breakpoint0
            unique_target1 = total_unique_targeted1/total_unique_breakpoint1

            unique_target_odds = unique_target0/unique_target1

            try:
                del_target_odds = \
                    (d0.total_del/d0.total_targeted_del_breakpoints)/(d1.total_del/d1.total_targeted_del_breakpoints)
            except ZeroDivisionError:
                del_target_odds = 0
            try:
                udel_target_odds = \
                    (d0.unique_targeted_del_breakpoints / d0.total_unique_del) / (d1.unique_targeted_del_breakpoints /
                                                                                  d1.total_unique_del)
            except ZeroDivisionError:
                udel_target_odds = 0
            try:
                ins_target_odds = \
                    (d0.total_targeted_ins_breakpoints/d0.total_ins)/(d1.total_targeted_ins_breakpoints/d1.total_ins)
            except ZeroDivisionError:
                ins_target_odds = 0
            try:
                uins_target_odds = \
                    (d0.unique_targeted_ins_breakpoints / d0.total_unique_ins) / (d1.unique_targeted_ins_breakpoints /
                                                                                  d1.total_unique_ins)
            except ZeroDivisionError:
                uins_target_odds = 0

            total_targeted_odds_ratio_list.append(total_target_odds)
            unique_targeted_odds_ratio_list.append(unique_target_odds)
            total_targeted_del_odds_ratio_list.append(del_target_odds)
            total_targeted_ins_odds_ratio_list.append(ins_target_odds)
            unique_targeted_del_odds_ratio_list.append(udel_target_odds)
            unique_targeted_ins_odds_ratio_list.append(uins_target_odds)

            odds_string += \
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}" \
                "\t{}\t{}\t{}\t{}\t{}\n"\
                .format(total_target_odds, unique_target_odds, del_target_odds, udel_target_odds, ins_target_odds,
                        uins_target_odds, total_breakpoint0, d0.total_del, d0.total_ins, total_targeted0,
                        d0.total_targeted_del_breakpoints, d0.total_targeted_ins_breakpoints, total_unique_breakpoint0,
                        d0.total_unique_del, d0.total_unique_ins, total_unique_targeted0,
                        d0.unique_targeted_del_breakpoints, d0.unique_targeted_ins_breakpoints, total_breakpoint1,
                        d1.total_del, d1.total_ins, total_targeted1, d1.total_targeted_del_breakpoints,
                        d1.total_targeted_ins_breakpoints, total_unique_breakpoint1, d1.total_unique_del,
                        d1.total_unique_ins, total_unique_targeted1, d1.unique_targeted_del_breakpoints,
                        d1.unique_targeted_ins_breakpoints)

        odds_labels = "Total Targeted\tUnique Targeted\tDel Targeted\tUnique Del Targeted\tIns Targeted\t" \
                      "Unique Ins Targeted\tSample_0 Total\tSample_0 tDel\tSample_0 tIns\tSample_0 Targeted\t" \
                      "Sample_0 tDel Targeted\tSample_0 tIns Targeted\tSample_0 Unique\tSample_0 uDel\tSample_0 uIns\t"\
                      "Sample_0 uTargeted\tSample_0 uDel Targeted\tSample_0 uIns Targeted\tSample_1 Total\t" \
                      "Sample_1 tDel\tSample_1 tIns\tSample_1 Targeted\tSample_1 tDel Targeted\t" \
                      "Sample_1 tIns Targeted\tSample_1 Unique\tSample_1 uDel Targeted\tSample_1 uIns Targeted\n"

        total_odds_mean = round(scipy.mean(total_targeted_odds_ratio_list), 2)
        del_odds_mean = round(scipy.mean(total_targeted_del_odds_ratio_list), 2)
        ins_odds_mean = round(scipy.mean(total_targeted_ins_odds_ratio_list), 2)

        unique_odds_mean = round(scipy.mean(unique_targeted_odds_ratio_list), 2)
        unique_del_odds_mean = round(scipy.mean(unique_targeted_del_odds_ratio_list), 2)
        unique_ins_odds_mean = round(scipy.mean(unique_targeted_ins_odds_ratio_list), 2)

        total975 = numpy.percentile(total_targeted_odds_ratio_list, 97.5, interpolation='linear')
        total25 = numpy.percentile(total_targeted_odds_ratio_list, 2.5, interpolation='linear')

        del975 = numpy.percentile(total_targeted_del_odds_ratio_list, 97.5, interpolation='linear')
        del25 = numpy.percentile(total_targeted_del_odds_ratio_list, 2.5, interpolation='linear')

        ins975 = numpy.percentile(total_targeted_ins_odds_ratio_list, 97.5, interpolation='linear')
        ins25 = numpy.percentile(total_targeted_ins_odds_ratio_list, 2.5, interpolation='linear')

        unique_total975 = numpy.percentile(unique_targeted_odds_ratio_list, 97.5, interpolation='linear')
        unique_total25 = numpy.percentile(unique_targeted_odds_ratio_list, 2.5, interpolation='linear')

        unique_del975 = numpy.percentile(unique_targeted_del_odds_ratio_list, 97.5, interpolation='linear')
        unique_del25 = numpy.percentile(unique_targeted_del_odds_ratio_list, 2.5, interpolation='linear')

        unique_ins975 = numpy.percentile(unique_targeted_ins_odds_ratio_list, 97.5, interpolation='linear')
        unique_ins25 = numpy.percentile(unique_targeted_ins_odds_ratio_list, 2.5, interpolation='linear')

        outstring = "Permutation Analysis Module v{}; {} Type Permutations run {}\n" \
                    "Target File:\t{}\nSegCopy File:\t{}\n\n" \
                    "\tTotalOddsMean\tUniqueOddsMean\tTotal 97.5\tTotal 2.5\tUnique 97.5\tUnique 2.5\n" \
                    "Total\t{}\t{}\t{}\t{}\t{}\t{}\nDel\t{}\t{}\t{}\t{}\t{}\t{}\nIns\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                    "\n\n{}\n{}" \
            .format(__version__, self.args.Permutation_Type, date.today().strftime("%Y-%m-%d"), self.args.Target_File,
                    self.args.Segment_File, total_odds_mean, unique_odds_mean, total975, total25, unique_total975,
                    unique_total25, del_odds_mean, unique_del_odds_mean, del975, del25, unique_del975, unique_del25,
                    ins_odds_mean, unique_ins_odds_mean, ins975, ins25, unique_ins975, unique_ins25, odds_labels,
                    odds_string)

        outfile = open("{0}{1}_odds_ratios.txt".format(self.args.Working_Folder, self.args.Job_Name), 'w')
        outfile.write(outstring)
        outfile.close()
        self.log.info("Sample Permutation Complete")

        return
        #
        # ratio_mean_list = []
        # ratio_std_list = []
        # ratio_list = []
        # odds_ratio_list = []
        # outstring = ""
        #
        # # Format data for output file.
        # for sub_group in natsort.natsorted(cell_permutation_data_dict):
        #     for key, values in cell_permutation_data_dict[sub_group].items():
        #         if key == "bp":
        #             break_point_mean = int(round(scipy.mean(values)))
        #             break_point_std = round(scipy.std(values), 2)
        #             break_point_median = int(round(scipy.median(values)))
        #         elif key == "intsect":
        #             intersect_mean = int(round(scipy.mean(values)))
        #             intersect_std = round(scipy.std(values), 2)
        #             intersect_median = int(round(scipy.median(values)))
        #         elif key == "bp/intsect":
        #             ratio_mean = scipy.mean(values)
        #             ratio_std = scipy.std(values)
        #             ratio_list.append(values)
        #
        #     ratio_mean_list.append(ratio_mean)
        #     ratio_std_list.append(ratio_std)
        #
        #     outstring += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}"\
        #         .format(break_point_mean, break_point_median, break_point_std, intersect_mean, intersect_median,
        #                 intersect_std)
        #     outstring += "\t"
        #
        # for l1, l2 in zip(ratio_list[0], ratio_list[1]):
        #     odds_ratio_list.append(l1/l2)
        #
        # t = stats.t.interval(0.95, df=self.freq_calc_iterations-1, loc=scipy.mean(odds_ratio_list),
        #                      scale=scipy.std(odds_ratio_list) / numpy.sqrt(self.freq_calc_iterations))
        #
        # pval = stats.ttest_1samp(odds_ratio_list, 1)
        #
        # outstring += "{0}\t{1}\t{2}\t{3}\t{4}\n"\
        #     .format(round(scipy.mean(odds_ratio_list), 2), round(scipy.std(odds_ratio_list), 2), round(t[0], 2),
        #             round(t[1], 2), pval[1])
        #
        # for v in odds_ratio_list:
        #     outstring += "{0}\n".format(v)
        #
        # outfile.write(outstring)
        # outfile.close()
        #
        # print("Permutation Analysis of Samples Complete.")
        #
        # return

    def bin_sizing(self):
        """
        This will re-bin the genome space to be consistent with the breakpoint definition or any other bin re-sampling
        desired.
        :return:
        """

        self.log.info("Begin Re-Binning the Genome Space.")
        new_list = []
        seg_num = 0

        for chrom in natsort.natsorted(self.seg_analyzer.chrom_list):
            self.log.debug("Binning Chromosome {0}".format(chrom))

            # Some chromosomes have no segments.
            try:
                chrom_slice = \
                    self.seg_analyzer.seg_copy_array[self.seg_analyzer.seg_copy_array[:, 1] == chrom.encode()]
                seg_count = chrom_slice.shape[0]
                coord_start = int(chrom_slice[0, 2])
            except IndexError:
                continue

            for i in range((seg_count-1)):
                if (i+1) < seg_count and (i+1) % int(self.args.Combine_Segments) == 0:
                    coord_stop = int(chrom_slice[i, 3])
                    new_list.append([seg_num, chrom.encode(), coord_start, coord_stop])

                    coord_start = int(chrom_slice[i+1, 2])
                    seg_num += 1

        self.log.info("Genome Space Successfully Re-Binned.")

        return numpy.array(new_list, dtype='object')

    def target_mapping(self):
        """
        This method maps the target coordinates onto the genome space coordinates.  Also prints a file of the map
        coordinates if the user so desires.  Returns a numpy array of the map.
        :return:
        """

        map_list = []
        self.bin_tracking_array = self.seg_analyzer.bin_tracking_array
        self.log.info("Spawning {0} jobs to begin building Target_Bed_Map_Array for permutation analysis."
                      .format(self.args.Spawn))

        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        for lst in p.starmap(self.sub_target_mapping,
                             zip(itertools.repeat(self.bin_tracking_array), itertools.repeat(self.target_bed_array),
                                 itertools.repeat(self.args), self.seg_analyzer.chrom_list)):

            map_list.extend(lst)

        map_list.sort(key=lambda x: x[0])

        if eval(self.args.Map_File):
            self.log.info("Writing Map File")
            file_data = ""
            map_file = open("{0}{1}_{2}_mapfile.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, self.args.Cell_Name), 'w')
            map_file.write("Chrom\tstart\tstop\trefBinID\ttargetBinID\ttargetCount\n")

            for row in sorted(map_list, key=itemgetter(0)):

                coord_start = int(self.bin_tracking_array[self.bin_tracking_array[:, 0] == row[0]][0, 2])
                coord_stop = int(self.bin_tracking_array[self.bin_tracking_array[:, 0] == row[0]][0, 3])
                chrom = self.bin_tracking_array[self.bin_tracking_array[:, 0] == row[0]][0, 1].decode()
                r_count = len(row[1])
                file_data += ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
                              .format(chrom, coord_start, coord_stop, row[0], row[1], r_count))

            map_file.write(file_data)
            map_file.close()
            self.log.info("Map File Written")

        self.log.info("Target_Bed_Map_Array built.")
        return numpy.array(map_list, dtype='object')

    @staticmethod
    def sub_target_mapping(bin_tracking_array, target_bed_array, args, chrom):
        """
        This is the subprocess that computes the mapping coordinates for a target bed file on a single chromosome.
        :param bin_tracking_array:
        :param target_bed_array:
        :param args:
        :param chrom:
        :return:
        """
        log = Tool_Box.Logger(args, parellel_id=chrom)
        log.debug("Traversing chromosome {0}".format(chrom))

        map_list = []

        seg_count = bin_tracking_array[bin_tracking_array[:, 1] == chrom.encode()].shape[0]
        target_row_count = target_bed_array[target_bed_array[:, 1] == chrom.encode()].shape[0]
        chrom_slice = bin_tracking_array[bin_tracking_array[:, 1] == chrom.encode()]
        target_slice = target_bed_array[target_bed_array[:, 1] == chrom.encode()]

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

        Tool_Box.delete(["{}{}_{}.log".format(args.Working_Folder, args.Job_Name, chrom)])

        return map_list

    def permute_data(self):
        """
        Permutations based on random selection of segments.  --Group_Size is derived from actual number of breakpoints
        for cell type.
        :return:
        """

        if not eval(self.args.Segment_Permutation_File):
            """
            User has selected no permutation outputs and should not be here so kick them out.
            """
            self.log.error("--Segment_Permutation_File must be set True")
            return

        cell_list = []

        for cell in self.args.Cell_Name.strip().split(','):
            cell_list.append(cell)

        """
        This block is designed to do Segment Permutation Analysis.  Randomly select groups of
        self.permutation_group_size from the Genome space, in this case the list of ID's from the
        self.bin_tracking_list, and intersect them with the Target Bed data.
        """
        self.log.info("Spawning {0} jobs to process {1} iterations each for segment permutation analysis."
                      .format(self.args.Spawn, self.args.Iteration_Count))

        self.target_bed_map_array = self.target_mapping()
        selection_space = self.bin_tracking_array[:, 0]
        intersect_space = self.target_bed_map_array[:, 0].tolist()

        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        p.starmap(self.intersection_iteration, zip(itertools.repeat(selection_space),
                                                   itertools.repeat(intersect_space),
                                                   itertools.repeat(self),
                                                   cell_list))

        self.log.info("Segment Permutation jobs done.  Compile any temporary data files into final results.")

    @staticmethod
    def intersection_iteration(selection_space, intersect_space, self, cell):
        """
        Controls the number of iterations and collects the data returned from the intersection function.

        :param selection_space:
        :param intersect_space:
        :param selection_space:
        :param self:
        :param cell:
        :return:
        """

        target_bed_map_array = self.target_bed_map_array

        iteration_count = 0
        tmp_intersection_list = []

        log = self.log

        # Get the number of unique breakpoints.
        self.seg_analyzer.break_points(cell=cell, permutation=True)
        unique_permutation_group_size = self.seg_analyzer.unique_permutation_group_size
        total_permutation_group_size = self.seg_analyzer.total_permutation_group_size

        # Retrieve a namedtuple of the intersect data
        total_targeted_data = [int(self.args.Total_Targeted), total_permutation_group_size-int(self.args.Total_Targeted)]
        unique_targeted_data = [int(self.args.Unique_Targeted), unique_permutation_group_size-int(self.args.Unique_Targeted)]

        log.info("Unique Combination Group Size {} for Segment Permutation: {}"
                 .format(cell, unique_permutation_group_size))
        log.info("Total Combination Group Size {} for Segment Permutation: {}"
                 .format(cell, total_permutation_group_size))

        permuted_sum_file = "{0}{1}_{2}_seg_random_bin_total.bed" \
            .format(self.args.Working_Folder, self.args.Job_Name, cell)

        permuted_sum_file_out = open(permuted_sum_file, "w")
        sum_file_data = "Total Targeted Segments\tTotal Odds Ratio\tTotal Fisher Ecact p_Value\t" \
                        "Unique Targeted Segments\tUnique Odds Ratio\tUnique Fisher Exact p_Value\n"

        while iteration_count < int(self.args.Iteration_Count):
            iteration_count += 1
            # numpy.random.shuffle(selection_space)
            u = numpy.intersect1d(numpy.random.choice(selection_space, unique_permutation_group_size, replace=False),
                                  intersect_space)
            t = numpy.intersect1d(numpy.random.choice(selection_space, total_permutation_group_size, replace=True),
                                  intersect_space)

            if len(u) > 0:
                for seg in u:
                    tmp_intersection_list.append(str(seg))

            t_intersect = len(t)
            u_intersect = len(u)
            t_oddsratio, t_pvalue = \
                scipy.stats.fisher_exact([total_targeted_data, [t_intersect, total_permutation_group_size-t_intersect]])

            u_oddsratio, u_pvalue = \
                scipy.stats.fisher_exact([unique_targeted_data, [u_intersect, unique_permutation_group_size-u_intersect]])

            sum_file_data += "{}\t{}\t{}\t{}\t{}\t{}\n"\
                .format(t_intersect, t_oddsratio, t_pvalue, u_intersect, u_oddsratio, u_pvalue)

            if iteration_count % int(self.args.Prog_Check) == 0:
                log.debug("{0} iterations of {1} for Segment Permutation Analysis using file {2}"
                          .format(iteration_count, int(self.args.Iteration_Count), permuted_sum_file))

        permuted_sum_file_out.write(sum_file_data)
        permuted_sum_file_out.close()

        log.debug("Processing {} data for output to segment permutation file.".format(cell))

        seg_permute_data = "For a given segment, how many times do targets intersect in {0} iterations?\n" \
                           "chrom\tstart\tstop\tSeg_ID\tTotal_Targets\tNumber_Per_Permutation\tbin_size\n"\
            .format(self.args.Iteration_Count)

        seg_counts = Counter(tmp_intersection_list)
        for seg in seg_counts:
            bed_tuple = \
                target_bed_map_array[numpy.where(target_bed_map_array[:, 0] == int(seg))][:, 1].tolist()[0]

            chrom = self.bin_tracking_array[self.bin_tracking_array[:, 0] == int(seg)][0, 1].decode()
            coord_start = self.bin_tracking_array[self.bin_tracking_array[:, 0] == int(seg)][0, 2]
            coord_stop = self.bin_tracking_array[self.bin_tracking_array[:, 0] == int(seg)][0, 3]
            counts = len(bed_tuple)*seg_counts[seg]
            count = counts / int(self.args.Iteration_Count)
            bin_size = int(coord_stop) - int(coord_start)
            seg_permute_data += ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
                                 .format(chrom, coord_start, coord_stop, seg, counts, count, bin_size))

        seg_permute_file = open("{0}{1}_{2}_seg_permute.txt"
                                .format(self.args.Working_Folder, self.args.Job_Name, cell), 'w')

        seg_permute_file.write(seg_permute_data)
        seg_permute_file.close()

        # Tool_Box.delete(["{}{}_{}.log".format(self.args.Working_Folder, self.args.Job_Name, cell)])
        log.info("{} Permutation Complete".format(cell))

        return
