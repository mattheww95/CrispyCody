#!/usr/bin/env python3

import pysam
import argparse
import sys
import os
import statistics
from collections import Counter
import time
import gc
import subprocess


class CodyAnalysis:
    """
    This prepare the coverage data from an input bam file,
    needs to calculate the mean median and mode
    """
    def __init__(self, bam_files, depth_threshold=5, genome_length=29903):
        self.bam_files = bam_files  # tuple of a bam and bai
        self.depth_threshold = depth_threshold
        self.genome_length = genome_length
        self.sample_cov = self.get_coverage_info()
        self.sample_mean = self.calculate_average()
        self.sample_median = self.calculate_median()
        self.breadth_of_coverage = self.breadth_coverage()
        self.sample_mode = None

    def get_coverage_info(self):
        """
        Select calculate coverage for each of the files selected.
        Need one sample per a call not all
        """
        bam = self.bam_files
        #bamfile = pysam.AlignmentFile(bam[0], "rb")
        # postion and depth as dict
        sample_cov_ = {}
        #for pileupcolumn in bamfile.pileup(contig="MN908947.3",
        #                                    filepath_index=bam[1],
        #                                   stepper='nofilter'):
        #    sample_cov_[int(pileupcolumn.pos)] = int(pileupcolumn.n)
        # TODO add samtools to installer trying subprocess call again instead
        #sample_cov = pysam.depth("-aa", f"{bam[0]}", catch_stdout=True)
        depth_cov = subprocess.Popen(["samtools", "depth", "-aa", str(bam[0])],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        depth_comm = depth_cov.communicate()[0].decode("utf-8")

        sample_cov = depth_comm.split("\n")
        for line in sample_cov:
            split_line = line.split("\t")
            if len(split_line) == 3:
                sample_cov_[int(split_line[1])] = int(split_line[2])
        sys.stdout.flush()
        # need to deal with empty positons
        # Genome is one based
        #for i in range(0, self.genome_length):
        #    if sample_cov_.get(i) is None:
        #        sample_cov_[i] = 0

        return sample_cov_


    def calculate_average(self):
        """
        Calculate the average depth of the genome, seperate from median as
        median needs to be sorted
        """
        avg_vals = [i for i in self.sample_cov.values()]
        avg_coverage = round(statistics.mean(avg_vals), 2)
        #sample_values = 0
        #for depth in self.sample_cov:
        #    if depth > self.depth_threshold:
        #        sample_values += depth
        #    else:
        #        pass

        #avg_coverage = sample_values / self.genome_length
        return avg_coverage

    def calculate_median(self):
        """
        Sort the data so that the median value can be calculated and the mode
        """
        # To sort the dicitionary by values
        # sorted_sample_dict = sorted(self.sample_cov.items(), key=lambda x:x[1], reverse=False)
        # sort dict witout position for now
        cov_vals = [i for i in self.sample_cov.values()]
        # sorted_sample_dict = sorted(cov_vals, reverse=False)
        # print(sorted_sample_dict) reutrns list of tupls
        # get list index for middle based on genomes length
        cov_median = statistics.median(cov_vals)
        #print(self.sample_cov)
        return cov_median

    def breadth_coverage(self):
        """
        Breadth of coverage is depth at each position over threshold, divided by length
        """
        vals_over_threshold = [i for i in self.sample_cov.values() if i >= self.depth_threshold]
        breadth_cov = round(len(vals_over_threshold) / float(self.genome_length) * 100, 2)
        return breadth_cov



class UserInteraction(object):
    """
    Process the command line arguments
    """

    def __init__(self):
        self.parameters = self.get_args()
        self.run_data()

    def get_args(self):
        parser = argparse.ArgumentParser(description="Calculate the breadth of"\
                                         " coverage")
        parser.add_argument('-b', '--bam', help="Specify an input bam file",
                            action='store', default=None)
        # not index arg this time, assumes it is there
        parser.add_argument('-t', '--threshold', help='Set the threshold for' \
                            'depth of pileup coverage default value is 5',
                            type=int, default=5, action='store')
        parser.add_argument('-d', '--directory_scan', action='store',
                            help="A path for directories to scan", default=None)
        parser.add_argument('-l', '--genome_length', help="Seth the length of the genome to be calculated, default is 29903", default=29903, type=int)

        args = parser.parse_args()
        return args

    @staticmethod
    def create_log(file_handle):
        """
        Create a log file for the data to run, takes a file handle
        """
        if os.path.isfile(file_handle):
            file_abs_path = os.path.abspath(file_handle)
            file_directory = os.path.dirname(file_abs_path)
            # log name
        else:
            file_directory = os.path.abspath(file_handle)

        log_path_out = os.path.join(file_directory, "Crispy_Cody.csv")
        return log_path_out


    def run_data(self):
        """
        This part of the program will batch out the actual analyses to perform
        based on cmd args passed
        """
        depth = self.parameters.threshold
        length = self.parameters.genome_length

        # decide if scanning a directory or not
        dir_scan_path = self.parameters.directory_scan
        single_bam = self.parameters.bam  # this was bam_files[0] before
        if dir_scan_path is None and single_bam is None:
            sys.exit("No files selected.")

        results = []
        if dir_scan_path is not None and os.path.exists(dir_scan_path):
            log_file = UserInteraction.create_log(dir_scan_path)
            directory_analysis = PrepareFiles(dir_scan_path)
            for analysis in directory_analysis.files:
                results.append(analysis)
        elif single_bam is not None and os.path.exists(single_bam):
            # Check for index, refactor into function later
            log_file = UserInteraction.create_log(single_bam)
            index_file = os.path.abspath(single_bam) + ".bai"
            if os.path.isfile(index_file):
                files_tupe = (single_bam, index_file)
                single_result = CodyAnalysis(files_tupe)
            else:
                base_file = os.path.abspath(single_bam)
                # create index
                indexing = subprocess.Popen(["samtools", "index", "-b", base_file])
                indexing.wait()
                single_result = CodyAnalysis(files_tupe)
            results.append(single_result)
        logfile = open(log_file, 'w')
        logfile.write("Bam_file,Percent_coverage,avg_depth,median_coverage\n")
        for result in results:
            logfile.write(f"{result.bam_files[0]},{result.sample_mean},"\
            f"{result.breadth_of_coverage},{result.sample_median}\n")
        logfile.close()
        print("Done")



class PrepareFiles:
    """
    Prepare the files for analysis
    Input will be the directory but need to initialize the class with the bams
    """

    def __init__(self, directory_bam):
        self.files = self.analyse_directory(directory_bam)

    def analyse_directory(self, directory_bam):
        """
        :param: directory_bams: The folder to find the bams in
        :return: Return a CodyAnalysis class in list
        """

        directory_bams = [os.path.join(directory_bam, i) for i in
                            os.listdir(directory_bam) if i[-3:] == "bam"]
        # check for no samples in the directory
        if len(directory_bams) == 0:
            sys.exit(f"No bam files found in {directory_bam}")

        # Find the bais
        files_to_return = []
        for bam in directory_bams:
            #index_file = os.path.join(bam, ".bai")
            index_file = os.path.abspath(bam) + ".bai"
            if os.path.isfile(index_file):
                file_tupe = (bam, index_file)
                files_to_return.append(file_tupe)
            else:
                base_file = os.path.abspath(bam)
                # basefile to string
                str_base = str(base_file)
                # Memory leak in pysam, buffers being over read trying a
                # garbage collection first, slows down performance hopefully solves
                # the issue though
                #gc.collect()
                #pysam.index(str_base, catch_stdout=False, catch_stderr=False)
                # Tired of pysam and its memory management
                # Switching to a subprocess
                # TODO add samtools to install
                indexing = subprocess.Popen(["samtools", "index", "-b", str_base])
                indexing.wait()
                # Adding sleep as hopefully only time to close the buffer is needed
                file_tupe = (bam, bam + ".bai")
                files_to_return.append(file_tupe)

        analysed_samples_return = []
        for sample_set in files_to_return:
            cur_obj = CodyAnalysis(sample_set)
            analysed_samples_return.append(cur_obj)

        return analysed_samples_return


def tests():
    data_tup = ("../data/test.bam", "../data/test.bam.bai")
    cody_obj2 = PrepareFiles("../data/create_index_test")
    print("End of index test, first")
    # All indexing must be done before any analysis or else seffault occurs
    cody_obj = CodyAnalysis(data_tup)
    print(cody_obj.sample_mean, cody_obj.sample_median,
          cody_obj.breadth_of_coverage)
    print("End of single sample Test")
    directory_test = PrepareFiles("../data/")
    for i in directory_test.files:
        print(i.sample_mean, i.sample_median, i.breadth_of_coverage)
    many_dir_test = PrepareFiles("../data/many_bams")
    for i in many_dir_test.files:
        print(i.sample_mean, flush=True)
    print("End of many bam test")
    print("End of directory test.")
    empty_dir_test = PrepareFiles("../data/empty_directory")


if __name__=="__main__":
    UserInteraction()
