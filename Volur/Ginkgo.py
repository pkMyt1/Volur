"""
Genome_Bin_Generator v0.5.0
    June 9, 2017
    Dennis A. Simpson
    Made changes to use some of the utilities from Tool_Box.  Added output that contains data about the coverage depth
    and breadth.  Modified the code so the variable bin search is faster now.
Genome_Bin_Generator v0.4.0
    February 13, 2017
    Dennis A. Simpson
    Have bin_mapper function to reverse engineer the bin data that Ginkgo uses.  Many tweaks on the Variable Bin
    Generator class for speed and efficiency.  Reduced the minimum bin size equal the step size.  Overall this is
    working quite well.  Still takes several hundred cpu hours to run.  Could benefit from Cythonizing.
Genome_Bin_Generator v0.2.2
    January 16, 2017
    Dennis A. Simpson
    Code uses multiprocessing now.  Found a bug where zero length chromosome names in the index file were causing a
    divide by zero error.  Fixed by adding the Try near line 104.
Genome_Bin_Generator v0.1.0
    January 12, 2017
    Dennis A. Simpson
    This will produce a file of custom bin sizes for Ginkgo.
    Copyright 2017
"""

import os
import pathlib
import pysam
import multiprocessing
import itertools
import collections
import natsort
import csv
import pyfaidx
import statistics
import Valkyries.Tool_Box as Tool_Box
import Valkyries.BAM_Tools as Bam_Tools

__version__ = "0.5.0"
__author__ = "Dennis A. Simpson"
__package__ = 'Mimir'


def bin_mapper(o):
    """
    This module is designed to reverse engineer the bins for Ginkgo starting with just bin coordinates.
    :param o:
    :return:
    """
    infile = list(csv.reader(open(o.input_file), delimiter='\t'))
    filepath = str(pathlib.PurePath(pathlib.PurePath(o.input_file).parent, ""))
    outfile_gc_content_name = "{0}{1}GC_variable_25000.{2}".format(filepath, os.sep, o.Aligner)

    outfile_bounds_name = "{0}{1}bounds_variable_25000.{2}".format(filepath, os.sep, o.Aligner)

    outfile_gc_content = open(outfile_gc_content_name, 'w')
    outfile_bounds = open(outfile_bounds_name, 'w')

    refseq_file = pyfaidx.Fasta(o.Ref_Seq)
    bounds_count = 0
    gc_content = ""
    bin_stop = 0
    prev_chrom = ""

    first_line = True
    for line in infile:
        if len(line) < 2:
            raise SystemExit("line {0} is too short".format(bounds_count))

        chrom = line[0]
        bounds_count += 1
        if first_line:
            bin_start = 1
            prev_chrom = chrom
            first_line =False
            print("Begin Processing {0}".format(chrom))
        else:
            bin_start = bin_stop + 1

        bin_stop = int(line[1])

        t2 = str(refseq_file.get_seq(name=chrom, start=bin_start, end=bin_stop))

        if prev_chrom == chrom:
            gc_content += "{0}\n".format((t2.count("G") + t2.count("C")) / (bin_stop - bin_start))
        else:
            outfile_bounds.write("{0}\t{1}\n".format(prev_chrom, bounds_count))
            outfile_gc_content.write("{0}".format(gc_content))
            gc_content = "{0}\n".format((t2.count("G") + t2.count("C")) / (bin_stop - bin_start))
            print("\t{0} Finished.\n Begin Processing {1}".format(prev_chrom, chrom))
            prev_chrom = chrom

    outfile_bounds.write("{0}\t{1}\n".format(prev_chrom, bounds_count))
    outfile_gc_content.write("{0}".format(gc_content))
    refseq_file.close()
    outfile_gc_content.close()
    outfile_gc_content.close()


class VariableBinGenerator:
    """
    This class generates custom bins for Ginkgo analysis.
    """
    def __init__(self, o):

        self.o = o
        self.algn_dict = Bam_Tools.BamTools(o)
        self.coverage_dict = {"chrom":"tuple"}
        self.module_director()

    def module_director(self):

        algn_list = list(self.algn_dict.build_dictionary(self.o.BAM_File).items())
        options_dictionary = self.o._asdict()
        p = multiprocessing.Pool(int(self.o.spawn))  # Define max number of parallel jobs.
        coverage = Tool_Box.CoverageCalculator(self.o.BAM_File)  # Initialize this object for use in parallel jobs.

        print("Spawning \033[1;35m{0}\033[m parallel jobs to process genomic bins." .format(self.o.spawn))

        # Determine starting read counts for each bin as well as read depth and breadth for each chromosome.
        data_bundle1 = options_dictionary, coverage
        tmp_dict_list = p.starmap(self.read_mapping, zip(itertools.repeat(data_bundle1), algn_list))

        for dict in tmp_dict_list:
            self.coverage_dict.update(dict)
        del tmp_dict_list

        print("-->Mapping complete. Begin Generating Custom Bins")
        # Calculate custom bins.
        data_bundle2 = self.coverage_dict, options_dictionary
        p.starmap(self.create_genomic_bins, zip(itertools.repeat(data_bundle2), algn_list))

        print("***All parallel jobs complete. Begin compiling data.***")

        # Initialize our output files.
        filepath = str(pathlib.PurePath(pathlib.PurePath(self.o.BAM_File).parent, ""))
        initial_bin_size = float(self.o.Bin_Size) * 1000

        outfile_name = "{0}{1}variable_{2}_{3}.{4}".format(
            filepath, os.sep, int(initial_bin_size), self.o.Read_Length, self.o.Aligner)
        outfile_name2 = "{0}{1}variable_{2}_{3}_{4}_stats.txt".format(
            filepath, os.sep, int(initial_bin_size), self.o.Read_Length, self.o.Aligner)
        outfile_gc_content_name = "{0}{1}GC_variable_{2}_{3}.{4}" .format(
            filepath, os.sep, int(initial_bin_size), self.o.Read_Length, self.o.Aligner)
        outfile_bounds_name = "{0}{1}bounds_variable_{2}_{3}.{4}".format(
            filepath, os.sep, int(initial_bin_size), self.o.Read_Length, self.o.Aligner)

        outfile_bin = open(outfile_name, 'w')
        outfile_stats = open(outfile_name2, 'w')
        outfile_gc_content = open(outfile_gc_content_name, 'w')
        outfile_bounds = open(outfile_bounds_name, 'w')

        outstring_bin = "CHR\tEND\n"
        bounds_count = 0
        outfile_bin.write(outstring_bin)
        tmp_file_list = []

        # Setup coverage output.
        outstring_stats = "CHR\tCoverage_Depth\tCoverage_Breadth"
        for chrom in natsort.natsorted(self.coverage_dict):
            outstring_stats += "\n{0}\t{1}\t{2}"\
                .format(chrom, self.coverage_dict[chrom][3], self.coverage_dict[chrom][4])
        outstring_stats += "\n\nCHR\tBinStart\tBinStop\tBinSize\tActualCounts\tEstimatedCounts\n"

        outfile_stats.write(outstring_stats)
        print("-->Output files initialized, begin reading temp data files.")

        for chrom in natsort.natsorted(Tool_Box.chromosomes(self.o.Species, self.o.chrY)):
            print("   -->Formatting output data for \033[1;35m{0}\033[m.".format(chrom))
            outstring_bin = ""
            outstring_stats = ""
            outstring_gc_content = ""

            infile_bin_name = "{0}{1}{2}_bin.tmp".format(filepath, os.sep, chrom)
            infile_stat_name = "{0}{1}{2}_stats.tmp".format(filepath, os.sep, chrom)

            if os.path.isfile(infile_bin_name) and os.path.isfile(infile_stat_name):
                tmp_bin_file = list(csv.reader(open(infile_bin_name), delimiter='\t'))
                tmp_stats_file = list(csv.reader(open(infile_stat_name), delimiter='\t'))

                tmp_file_list.extend([infile_bin_name, infile_stat_name])

                for bin_file_line, stats_file_line in zip(tmp_bin_file, tmp_stats_file):
                    bounds_count += 1
                    outstring_bin += "{0}\t{1}\n".format(bin_file_line[0], bin_file_line[1])
                    outstring_stats += "{0}\n".format("\t".join(stats_file_line))
                    outstring_gc_content += "{0}\n".format(bin_file_line[2])

                outfile_bin.write(outstring_bin)
                outfile_stats.write(outstring_stats)
                outfile_gc_content.write(outstring_gc_content)
                outfile_bounds.write("{0}\t{1}\n".format(chrom, bounds_count))

                print("    -->\033[1;35m{0}\033[m complete, data written to output files.".format(chrom))
            else:
                print("\033[1;31mCAUTION:\033[m\n\t-->{0} and/or {1} not found.  This may not be an error if the "
                      "number of mapped alignments was low".format(infile_bin_name, infile_stat_name))

        print("***Data compiled, begin disk cleanup.***")
        outfile_bin.close()
        outfile_stats.close()
        outfile_gc_content.close()
        outfile_bounds.close()
        Tool_Box.delete(tmp_file_list)

        return

    @staticmethod
    def read_mapping(data_bundle1, region):
        """
        This gives us the distribution of reads accross fixed bins.  This data is used to determine the target read
        counts in each variable bin.

        :param region:
        :param data_bundle1:
        :return:
        """

        options_dict = data_bundle1[0]
        c = data_bundle1[1]

        o = collections.namedtuple('options_file', options_dict.keys())(**options_dict)
        chrom = region[0]
        coverage_depth, coverage_breadth = c.coverage(region)

        bin_size = float(o.Bin_Size) * 1000
        bin_start = 1
        bin_stop = bin_size
        coverage_list = []
        coverage_dict = {}
        eof = False
        mark = False
        bam_file = pysam.AlignmentFile(o.BAM_File)

        print("    -->Calculating statistics for bins on \033[1;35m{0}\033[m.".format(chrom))
        while not eof:
            try:
                read_count = bam_file.count(reference=chrom, start=bin_start, end=bin_stop)
                if mark:
                    eof = True
            except ValueError:
                print("\n\033[1;31mWARNING:\033[m\n\tError in BAM file processing {0}, start: {1}, stop: {2}\n"
                      .format(chrom, bin_start, bin_stop))
                break

            coverage_list.append(read_count)

            # This block prevents running past the end of the chromosome and prevents errors when the bin size is
            # greater than the chromosome size.
            if bin_stop+bin_size > int(region[1][0]):
                if bin_size <= int(region[1][0]):
                    bin_start = int(region[1][0]) - bin_stop + 1
                else:
                    bin_start = 1

                bin_stop = int(region[1][0])
                mark = True
            else:
                bin_start = bin_stop+1
                bin_stop += bin_size

        coverage_dict[chrom] = statistics.mean(coverage_list), statistics.median(coverage_list), \
                               statistics.stdev(coverage_list), coverage_depth, coverage_breadth
        print("     -->Statistics calculated for bins on \033[1;35m{0}\033[m.".format(chrom))

        return coverage_dict

    @staticmethod
    def create_genomic_bins(data_bundle2, region):
        """
        This function is called by the multiprocessor.  It will process one chromosome at a time setting the size of
        the bins, stats about the bins, and GC content of the bins.

        :param data_bundle2:
        :param region:
        :return:
        """
        coverage_dict = data_bundle2[0]
        options_dict = data_bundle2[1]
        chrom = region[0]
        print("   --> Binning \033[1;35m{0}\033[m.".format(chrom))
        o = collections.namedtuple('options_file', options_dict.keys())(**options_dict)
        outstring_bin = ""
        outstring_stats = ""
        bins_found = 1
        kill_count = 0.0
        eof = False
        region_set = False
        step_count = 0

        initial_bin_size = float(o.Bin_Size) * 1000
        back_step = 0.1

        initial_bin_counts = int(coverage_dict[chrom][1])

        if initial_bin_counts > 1:
            initial_bin_counts = int(initial_bin_counts)
        else:
            initial_bin_counts = 1

        if initial_bin_counts < 10:
            print("\n\033[1;31mWARNING:\n\tEstimated read count per bin for {0} is <10.  "
                  "\n\tNot enough alignments to bin chromosome.\n\033[m".format(chrom))
            return
        elif initial_bin_counts < 100:
            print("\n\033[1;31mWARNING:\n\tEstimated read count per bin for {0} is {1}."
                  "\n\tWhen count is <100, quality of output in question.\n\033[m".format(chrom, initial_bin_counts))

        # Open files for reading.
        bam_file = pysam.AlignmentFile(o.BAM_File)
        refseq_file = pyfaidx.Fasta(o.Ref_Seq)

        bin_start = 1
        bin_stop = int(initial_bin_size)
        step = int(o.Step_Size)

        while not eof:
            try:
                t = bam_file.count(reference=chrom, start=bin_start, end=bin_stop)
            except OSError:
                print("\n\033[1;31mWARNING:\033[m\n\tError in BAM file processing {0}, start: {1}, stop: {2}\n"
                      .format(chrom, bin_start, bin_stop))
                break
            if (0.9 * initial_bin_counts) <= t <= (1.1 * initial_bin_counts):
                region_set = True

            # If we are under on the counts increase the bin size.
            elif t < (0.9 * initial_bin_counts):
                bin_stop += step
                step_count += 1

                # Going past the end of a BAM file region is not a good idea
                if bin_stop > region[1][0]:
                    bin_stop = int(region[1][0])-1
                    region_set = True
                    eof = True

            # If we are over on the counts then walk the bin size back.
            elif t > (1.1 * initial_bin_counts):

                bin_stop -= int(step * back_step)
                kill_count += 1

                # This keeps us from getting stuck in a loop between too small and too large.
                if bin_stop <= bin_start or kill_count >= 1/back_step:
                    kill_count = 0
                    bin_stop += int(step * back_step)
                    region_set = True

            # If the counts are just right, add the data to the output strings.
            if region_set:
                region_set = False
                bins_found += 1
                step_count = 0
                kill_count = 0
                t2 = str(refseq_file.get_seq(name=chrom, start=bin_start, end=bin_stop))
                gc_content = "{0}".format((t2.count("G") + t2.count("C")) / (bin_stop - bin_start))
                outstring_bin += "{0}\t{1}\t{2}\n".format(chrom, int(bin_stop), gc_content)
                outstring_stats += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    chrom, bin_start, int(bin_stop), int(bin_stop - bin_start), t, initial_bin_counts)

                bin_start = bin_stop + 1
                bin_stop += int(initial_bin_size)

                # Going past the end of a BAM file region is not a good idea
                if bin_stop > region[1][0]:
                    bin_stop = int(region[1][0]) - 1
                    if bin_start >= region[1][0]-1:
                        eof = True

                if bins_found % int(o.prog_check) == 0:
                    print("     -->Found \033[1;35m{0}\033[m bins of estimated \033[1;35m{1}\033[m for chromosome "
                          "\033[1;35m{2}\033[m."
                        .format(bins_found, int(region[1][0] / initial_bin_size), chrom))

        # Initialize and write data to our output files.
        filepath = str(pathlib.PurePath(pathlib.PurePath(o.BAM_File).parent, ""))
        outfile_bin_name = "{0}{1}{2}_bin.tmp".format(filepath, os.sep, chrom)
        outfile_stat_name = "{0}{1}{2}_stats.tmp".format(filepath, os.sep, chrom)

        outfile_bin = open(outfile_bin_name, "w")
        outfile_stats = open(outfile_stat_name, "w")

        outfile_bin.write(outstring_bin)
        outfile_stats.write(outstring_stats)

        outfile_bin.close()
        outfile_stats.close()
        bam_file.close()
        refseq_file.close()

        return

    # @staticmethod
    # def __build_dictionary(input_file):
    #     """
    #     Get the number of reads aligned for each chromosome and put this into a dictionary.  The "*" in position 0 are
    #     the total unmapped reads.  The idxstats format is region (Chr usually), region length, mapped, unmapped.
    #     :return:
    #     """
    #
    #     if len(input_file) < 3:
    #         raise SystemExit('\n\033[1;31mERROR:\033[m\tBAM file parameter missing from options file.')
    #     elif not os.path.isfile(input_file):
    #         raise SystemExit('\n\033[1;31mERROR:\033[m\tBAM file {0} not found'.format(input_file))
    #
    #     algn_dict = {}
    #
    #     for item in pysam.idxstats(input_file, split_lines=True):
    #         if not item.split("\t")[0] == "*":
    #             algn_dict[item.split("\t")[0]] = (int(item.split("\t")[1]), int(item.split("\t")[2]))
    #
    #     return algn_dict
