"""
Alignment_Launcher v0.4.1
    Dennis A. Simpson
    August 22, 2018
    Rename Ref_Seq variable to Aligner_Ref_Seq.  The module requires a reference sequence that is indexed to the
    specific aligner.  There is no reason this cannot all live in the same directory but UNC configures it differently.
Alignment_Launcher v0.4.0
    Dennis A. Simpson
    January 8, 2018
    Semantic versioning, pathos module for multiprocessor calls, testing mutithreading option for BWA.
Alignment_Launcher v0.75
    Dennis A. Simpson
    January 14, 2016
    This is in the process of becoming a generic module that can launch multiple types of aligners in parallel.
Alignment_Launcher v0.1
    Dennis Simpson
    August 13, 2015
    This module is designed to count the number of lines in the FASTQ files, break the file into a number of temp files
    equal to the --spawn parameter.  Each of these will be aligned separately then combined into a single SAM or BAM
    file.  Eventually I want this module to accept any aligner, the parameters for that aligner, and produce a sorted
    BAM file and the corresponding BAI file.
"""

import itertools
import ntpath
import pathos
import subprocess
from Valkyries import Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = "0.4.2"


class AlignmentLauncher:

    def __init__(self, args, log, paired_end):

        self.args = args
        self.paired_end = paired_end
        self.log = log

    def module_director(self, splitter_data):
        """
        As the name implies this function determines the size of the FASTQ files and coordinates writing and aligning
        them.
        :return:
        """

        p = pathos.multiprocessing.Pool(int(self.args.Spawn))
        p.starmap(self.run_aligner, zip(itertools.repeat(self.args),
                                        itertools.repeat(self.paired_end), splitter_data.fastq_file_list))

        return

    def run_aligner(self, fastq1, fastq2, bam_file=None):
        """
        This will run our aligner of choice.  Final output is a BAM file.
        :param fastq1:
        :param fastq2:
        :param bam_file:
        :return:
        """

        fq1_name = ntpath.basename(fastq1)
        fq2_name = None
        if self.paired_end:
            fq2_name = ntpath.basename(fastq2)

        sam_file = "{}{}.sam".format(self.args.Working_Folder, self.args.Job_Name)
        if bam_file is None:
            bam_file = fastq1.replace(".fastq", ".bam")

        if int(self.args.Spawn) > 2:
            threads = int(self.args.Spawn)
        else:
            threads = 1
        aligner_options = getattr(self.args, "Aligner_Options", "")
        if self.args.Aligner == "BWA":
            if self.args.BWA_Method == "mem":
                self.log.info("Begin alignment of {0} and {1} with BWA mem.".format(fq1_name, fq2_name))
                cmd = "bwa mem -t {0} {1} {2} {3} {4} > {5}" \
                    .format(threads, aligner_options, self.args.Aligner_Ref_Seq, fastq1, fastq2, sam_file)
                self.log.debug(cmd)
                subprocess.run([cmd], shell=True)

                self.log.info("Alignment complete, begin SAM to BAM conversion.")
                cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
                subprocess.run([cmd], shell=True)
                self.log.debug(cmd)

            elif self.args.BWA_Method == "aln":
                self.log.info("Begin alignment of {0} and {1} with BWA aln.".format(fq1_name, fq2_name))
                sai_file1 = fastq1.replace(".fastq", ".sai")
                sai_file2 = fastq2.replace(".fastq", ".sai")

                cmd1 = "bwa aln {0} {1} {2} {3} > {4}" \
                    .format(threads, aligner_options, self.args.Aligner_Ref_Seq, fastq1, sai_file1)
                cmd2 = "bwa aln {0} {1} {2} {3} > {4}" \
                    .format(threads, aligner_options, self.args.Aligner_Ref_Seq, fastq2, sai_file2)

                subprocess.run([cmd1], shell=True)
                subprocess.run([cmd2], shell=True)
                self.log.info("Alignment complete, begin SAI to SAM conversion.")

                if eval(self.paired_end):
                    cmd3 = "bwa sampe {0} {1} {2} {3} {4} > {5}" \
                        .format(self.args.Aligner_Ref_Seq, sai_file1, sai_file2, fastq1, fastq2, sam_file)

                elif not eval(self.paired_end):
                    cmd3 = "bwa {0} samse {1} {2} > {3}".format(self.args.Aligner_Ref_Seq, sai_file1, fastq1, sam_file)

                subprocess.run([cmd3], shell=True)

                self.log.info("SAI to SAM conversion complete. Begin SAM to BAM conversion.")
                cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
                subprocess.run([cmd], shell=True)

            else:
                self.log.error("No valid BWA alignment method provided")
                raise SystemExit(1)

        elif self.args.Aligner == "Bowtie2":
            if self.paired_end:
                sub_cmd = " -1 {0} -2 {1}" .format(fastq1, fastq2)
                message = "file {0} and {1}" .format(fq1_name, fq2_name)
            else:
                sub_cmd = " -U {0}" .format(fastq1)
                message = "file {0}" .format(fq1_name)

            if eval(self.args.local):
                bowtie2_local = " --local --ma {0}" .format(self.args.ma)
            else:
                bowtie2_local = ''

                self.log.info("Begin alignment of \033[1;35m{0}\033[m with Bowtie2" .format(message))

            # Construct command for aligner and call subprocess to run.
            cmd = "bowtie2 {0} --trim5 {1} --trim3 {2} -x {3} {4} -S {5}" \
                .format(bowtie2_local, self.args.trim5, self.args.trim3, self.args.Aligner_Ref_Seq, sub_cmd, sam_file)
            subprocess.run([cmd], shell=True)

            self.log.info("Alignment complete. Converting SAM format to BAM format")
            # Construct command for samtools and call subprocess to run.
            cmd = "samtools view -bh {0} -o {1}" .format(sam_file, bam_file)
            subprocess.run([cmd], shell=True)

        else:
            self.log.error("\033[1;31m***Warning:\033[m {0} not allowed.  Currently only BWA or Bowtie2 supported."
                           .format(self.args.Aligner))
            raise SystemExit(1)

        self.log.debug("File conversion process complete.  Removing temporary files.")
        Tool_Box.delete([fastq1, fastq2, sam_file, "{0}{1}_{2}.log"
                         .format(self.args.Working_Folder, self.args.Job_Name, fq1_name)])

        return True
