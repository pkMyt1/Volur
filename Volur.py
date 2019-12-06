#!usr/bin/env python3
"""
This handles Single Cell Sequence Analysis.
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2018
"""

import argparse
from argparse import RawTextHelpFormatter
import time
import sys
import os
import pathos
import Valkyries.Alignment_Launcher as Alignment_Launcher
import Valkyries.Version_Dependencies as VersionDependencies
import Volur.Permutation_Analysis as Permutation_Analysis
import Volur.Segment_Analyzer as Segment_Analyzer
import Valkyries.BAM_Tools as BAM_Tools
import Valkyries.Tool_Box as Tool_Box


__author__ = 'Dennis A. Simpson'
__version__ = '0.1.0'
__package__ = 'Volur'


def main(command_line_args=None):
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv

    parser = argparse.ArgumentParser(description="A package to process Single Cell Data.\n {0} v{1}"
                                     .format(__package__, __version__), formatter_class=RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()
    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)

    start_time = time.time()
    module_name = ""

    if eval(args.Permutation_Analysis):
        log.info("{0} v{1}; Module: Permutation_Analysis v{2} Beginning"
                 .format(__package__, __version__, Permutation_Analysis.__version__))
        options_parser.add_argument("--Target_Intersect_File", dest="Target_Intersect_File", default="False")
        options_parser.add_argument("--Breakpoint_Coordinate_File", dest="Breakpoint_Coordinate_File", default="False")
        args = options_parser.parse_args()
        pa = Permutation_Analysis.PermutationAnalysis(log, args)
        module_name = "Permutation_Analysis"

        if args.Permutation_Type == "Segment":
            pa.permute_data()
        elif args.Permutation_Type == "Sample":
            pa.cell_permutation()
        elif args.Permutation_Type == "Shuffle":
            pa.shuffle_analysis()
        else:
            log.error("No Permutation Analysis Type Selected")
            return

    elif eval(args.Single_Cell):
        index_list = Tool_Box.FileParser.indices(log, args.Index_File)
        paired_end = False
        fastq2 = getattr(args, "FASTQ2", None)
        if fastq2 is not None:
            paired_end = True

        # Now get to work.
        aligner = Alignment_Launcher.AlignmentLauncher(args, log, paired_end)
        # Aligner using multiple threads option
        bam_file = "{}{}.bam".format(args.Working_Folder, args.Job_Name)
        aligner.run_aligner(args.FASTQ1, fastq2, bam_file)

        log.info("Sending BAM file to Bam_Tools.")
        sorted_bamfile_name = BAM_Tools.sort(bam_file, int(args.Spawn), args.Compression_Level)
        BAM_Tools.index(sorted_bamfile_name)
        bamtool = BAM_Tools.BamTools(args, log, index_list)
        bamtool.demultiplex()
        bedfile_list = bamtool.bam_to_bed()

        # Get the BED files compressed.
        p = pathos.multiprocessing.Pool(int(args.spawn))
        p.starmap(Tool_Box.compress_files, zip(bedfile_list))

        module_name = "Single_Cell"

    elif eval(args.Segment_Analyzer):
        log.info("{0} v{1}; Module: Segment_Analyzer v{2} Beginning"
                 .format(__package__, __version__, Segment_Analyzer.__version__))

        output_validation = [eval(args.Breakpoint_Dist_File), eval(args.Breakpoint_Chrom_Dist_File),
                             eval(args.Chrom_Ploidy_File), eval(args.Aberration_Size_File),
                             eval(args.Target_Intersect_File), eval(args.Ensembl_Gene_File),
                             eval(args.Breakpoint_Coordinate_File)
                             ]

        if True not in output_validation:
            log.error("Segment Analyzer:  No output file type selected.")
            raise SystemExit(1)

        # Runs the analysis one cell type at a time sequentially.
        # This really is not the best approach.  A multiprocessor approach would likely be more efficient.
        cell_names = args.Cell_Name
        for cell in cell_names.strip().split(','):
            if eval(args.Target_Intersect_File) or eval(args.Breakpoint_Coordinate_File):
                if not os.path.isfile(args.Target_File):
                    log.error("{0} not found.  Check file name and path in Options File and try again."
                              .format(args.Target_File))
                    raise SystemExit(1)

            args.Cell_Name = cell

            # The way it is written this class must be initialized for each cell type.
            sa = Segment_Analyzer.SegmentAnalyzer(log, args)
            sa.chromosome_ploidy()
            sa.break_points()

            if eval(args.Aberration_Size_File):
                sa.size_output()

            if eval(args.Ensembl_Gene_File):
                # This processes the recurrent gene intersects by region.
                sa.genes()

            if eval(args.Target_Intersect_File) or eval(args.Breakpoint_Coordinate_File):
                sa.target_intersection()

        module_name = "Segment Analyzer"

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Volur {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


if __name__ == '__main__':
    main()
