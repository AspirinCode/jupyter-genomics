""" Module to test exec_megahit member methods. This module contains unit tests for exec_megahit.py. """

import os
import unittest
import logging
import shutil

from src import exec_megahit

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"

# input file contents
_TRIM_R1_FASTQ_STR = """@K00180:79:H3NG2BBXX:1:1115:6888:46135 1:N:0:GCCAAT+AGATCTCG
CGCTTCACCGCAATTGTGACCACCACCACAGCAGCAACCGGCAGCTTCGTTCTGGATCTTTTCGGCAGCCACTTCGTTGAAGTTCGGGAACTGGTTGGTTCCTAACAGGATTTCACGACGGGTAGCGACTGCTTTCTTGCGGGCTTCGTT
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<JJFJJJJJJJJJJJJJJJJJJJJJJJJJJJAAJFJJJJJJJJJJFJJJJJJJJJJJJJFJFFJJFJJJFJAFJJJJJJJJJJJJJJJ<J<<A-FFJJJ-7-AFJ
@K00180:79:H3NG2BBXX:1:2110:21298:16348 1:N:0:GCCAAT+AGATCTCG
ATAATTTGTATTTAAATATTATGTAGTGTAACACACAATTCTCATAAGATGTTGAAATTGTAGTATGTATTGATTCATAAGTCAGTTTTAAACTATTTCCTTACTTTCCTTCAATTGAAAGTATGATCATATATGTATTATCGTGTACAA
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAAAFFJ7<<FJJJJJ<FAAJA--AFJFFFF
@K00180:79:H3NG2BBXX:2:2222:8623:23751 1:N:0:GCCAAT+AGATCTCG
AAATGGGCTTATACGCGCTACAGAAGCGATAAAGGAGGTACCATGTGAACGCAGAAAAAACGAATATCCCAAACAGCGTGGATGCAACCCAAATTCCGGAATACGTCTTTGAATCGCTGGCACGGAGCCTGCTTCCGCTAACTCAGAAGT
+
AAFFFJJJJJJJJJJJJJJFJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJ<JJJJJJJJJJJJJJJFFJJ-AAJJJJJJJJJJJJAJJJJJJJJJJJ7FJJJJAFFFAJFJFJJJJJFFJJFA)<FJ77FJJJJJJJ)F-77<FFF<<
@K00180:79:H3NG2BBXX:2:1101:28970:1086 1:N:0:NCCAAT+NGATCTCG
AGGTCGAGCAGGATGCCTGCCAGCACCTTCTTGTAACGGATCGTATCGCCGCGCAGCGCGTCGGTAGTGACGTTGTCGATGTTCTTGATGAACACCACGTCGGCGTCGATCTCGGTCAGGTGTTCGATCAGCGCGCCGTGGCCCGCAGG
+
AAA<FJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJ<JF<JA<J<J-FFAJJJJ7FAJJJJJ<AJJ-FAJFJ7FAA<AFFAFFAJFFJ<FJJFJJJAAJF77F<-A-F-A-AF)A<AJ<-)7<<AA7A<<A<<7<F)7-7<<F-A77F
@K00180:79:H3NG2BBXX:2:2225:10247:10950 1:N:0:GCCAAT+AGATCTCG
GCCATAGCCGGGATATCCTGATACCGCATATCGCCGACTGGCAGCCGTCTGGCTTCCATCCGCTCATCCCAAAAAGGCGTCCCTTCTTGGCTGCGGCAGCATACGTAGAACCTATGTCGTCGTCGGAACAAAGTCCGTTCCACTCCGC
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJ<JFA-7AJJJJJFJJFJJJJA-AJ7JJJJFAFFJJJJJJFFJJJFJJJJJAFFJFJJJ7F<-AA7<)F<JAFJJJAFAJ77-F-A7AA7F"""
_TRIM_R2_FASTQ_STR = """@K00180:79:H3NG2BBXX:1:1115:6888:46135 2:N:0:GCCAAT+AGATCTCG
GGAAGCCATGTCTGCTGCTCTTGCCGGTGTCGATTCGATCACGGTTCGTCCGTTCGACAAGACTTATCAGACTCCGGACGACTTCTCCGAACGTATAGCCCGCAACCAG
+
AA<FFAJJJJJAJFFJJJJFFJ<FJJFJJJFFJJJJJJJJJJJJ<JJJJ7F<J-AJJJJJAFFJJJJJJJFJJJJJJJJJ-AJAFJJA<)<F<FFJ-7<FJJJFFJ)-<
@K00180:79:H3NG2BBXX:1:2110:21298:16348 2:N:0:GCCAAT+AGATCTCG
ATTTATATCTCTTTCCTAGAGGGTACAGTCTATGTCCTAGTCATCTATGTATTTTTCTTTTCACCTAGCACAATATCAGTTGCTAGTAGAATGAATTAAGAAACTGCTCCATACATTTTTTTTCTAATTTTCTGCTTTATCTAATGTTTC
+
AAAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJAJJ7FFFJJFJJJJFJJJJJJJJJAFJJJJJJJJJAJFJJJF7F7)AFFJFJFJJFJJJFJ77FFJF-AF-<FJ<F-A--<-AFFJ
@K00180:79:H3NG2BBXX:2:2222:8623:23751 2:N:0:GCCAAT+AGATCTCG
GCATTAGCACAGTTGGTAGCGCGCAACGTTCGCAATGTTGAGGTCAGGGGTTCGACTCCCCTATGCTCCACCAAAGCAAGATAATCCGAACCTTAGACCAATCGG
+
AAFFFJJFJJJJJJJJJJJJJJJJJJFJJJJJJJFAFJFJFJJAF<AJ<FF--7FAJJFFFJ7<<FAFJAJJJ-FFFA-FFJJ<-7FA7FAFA7F7JJFF--A<A
@K00180:79:H3NG2BBXX:2:1101:28970:1086 2:N:0:NCCAAT+NGATCTCG
CGTCACCTCCCACGCCTATCCCGAAGGCGCGCGCAAGGCCGTCGAGGAACATCTGGTCGAAGGGGCGGTCTATGCGGCCGCGAACGGC
+
<<<AA<J<-7JJJAJA7F<AFJ7A<777-77J-7-<-FF-7-7-777-AFA7AA--7AFF77--7A-77A-A-77-7AJAJA<AJF)<
@K00180:79:H3NG2BBXX:2:2225:10247:10950 2:N:0:GCCAAT+AGATCTCG
AGTTCCTGCGCCCGGAGGACCGGGTGCTGCTCATCGACGATTTCCTGGCCAACGGCAGCGCCCTACAGGGGCTTATCAAGCTGGCGGAGG
+
AAAFFJJJJJJJJJJJJJJJJJJJ-FFJJJAJJFJJJJJJFFJJJJJFA7A7FJFJ7JJJJJFF<JF<JF<F7--77<-AFFJJ)-AJF7"""
_TRIM_U_FASTQ_STR = """@K00180:79:H3NG2BBXX:2:2116:24322:14414 1:N:0:GCCAAT+AGATCTCG
GTTCGTTTGCCGTAGAAATCGCACGAAAAATCGGATTGCCGGAAGATGTCATTGCAGATGCTTCGGAAATCGTAGGAAGTGAATATATCAATGCCGACAAATACCTTCAGGACATCGTGCGCGACAAGCGTTACTGGGAAGGAAAACGTC
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJF<FAJJJJJJJJJJJJJJJJFJJJFFJJFFJJJFJJ<FFFFFJFAJJFJJ<FJJ<<JFJJJJFJAFJJ7-<A<<
@K00180:79:H3NG2BBXX:1:2127:15666:24982 1:N:0:GCCAAT+AGATCTCG
ACTGACGTGCGTCCTTGGGGCGCTGAAATTTTCACCCGCGATTTGGCTATGGTTCTCGGCATCAGTCTTGAGGAAGCCGAGGAATTAAAGCGCAATTCGGGGGAGTGCCGTCTGTCGCAGGTTATTGCCGGTGAAGTTGTACAGTTTGAA
+
AAAAFJFJJJJJJJJ7AJJJJJJJJJJJJJJJJJJJJJJJJFJJFFFJ<AJJFFJJJJFJJJJJJAJJJJJJJJJJJFJJJJJJJ77FFAFJJAJJFJFJ7FFJF-A<F7<J<JAAJF<AJAJJJJJFJJJ<A)<F-<<A<<7FAAJJ<J
@K00180:79:H3NG2BBXX:2:2208:1326:34037 1:N:0:GCCAAT+AGATCTCG
CAGAGGCAAGGCCTGAACCGCGAGGTTCCGACTGAAGGAAGCCCGAGGCAAACCATTGACCTTACGAACAGAAACGGTCTCAGGCATATGCTGTCGGGTAAGACTGCAAGACAAGTTGAAGCCCAAAGGCCACACGGAGATGGCAGTGTA
+
AAFFFJJJJJJJJJJJFJJJJJJJJJJF<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJAFFFJJJJJJJJJJJF<FFJJ7<AFJJJJJJJJJFJJJJJJJJJJJJF)<<FAAAJJJ7A)-<A7A77A-A-7
@K00180:79:H3NG2BBXX:2:1101:19116:1103 1:N:0:NCCAAT+NGATCTCG
CTGCTGTACCTGCTGGCTCTCCATTTTTTACCCAGCCCATCCAGCCATAGGTCTGTGCATGTACACGATAATATACATCATATAAATCTGCATTATTACCAGTCAGACGAATCTGAATTGCCTCCAGACGTTTTGACTGTCCACTTGTG
+
AAFFJJJJJJJJJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJFJFF7FJAF7FJJJJFJF-7F<FAAJJJJJJJ-<FFFFJJJJJJJAFFFF
@K00180:79:H3NG2BBXX:2:2116:9151:42636 1:N:0:GCCAAT+AGATCTCG
ATGATGATGAACACCACGATGGCCCATGCCGCCTCGGTGTTCAGGTCTGCCGTGGGGCTCCACAGGCCCACCACGCTGATGAGGTTGCAGCCGATGCTGAGTGCAAAGATAGCACCCACCAGCGGGATG
+
<AAFFJJJJJJJJ<JFJJJJJJJJJJJ<AFAJJJFJFFJJ-<JJJFFFFJJJJF<FJJJFJJJJJFJJJFJJF-AJJJ<-FFJJF-A7-FJAJ7FF<-AAJAF--FF<-77<-<<)<FFA7AAAAF))7"""


_MEGAHIT_PATH = '/usr/local/megahit/megahit'


class test_exect_trim(unittest.TestCase):
    """ Unit test exec_megahit methods """

    def setup_path(self):
        """ create strings corresponding to the temporary directories and files to be used in the unit tests """

        # trim input names
        _TRIM_R1_FQ_STR = 'trim_R1_20.fastq'
        _TRIM_R2_FQ_STR = 'trim_R2_20.fastq'
        _TRIM_U_FQ_STR = 'trim_U_20.fastq'

        # input file names
        _SAMPLE_R1_FQ_STR = 'readsInterleaved1.fa.gz'
        _SAMPLE_R2_FQ_STR = 'readsInterleaved2.fa.bz2'
        _SAMPLE_U_FQ_STR = 'readsInterleaved3.fa'

        # output file names
        _INTERMED_CONTIGS = 'intermediate_contigs'
        _FINAL_CONTIG = 'final.contigs.fa'
        _MEGAHIT_OPTS = 'opts.txt'

        # output subdirectories
        _TMP = 'tmp'
        _K21 = 'k21'

        # log file for the pipeline
        _PIPELINE_LOG = 'pipeline.log'

        # temporary directories to be used
        _UNITTEST_DIR_STR = 'megahit_unittest_temp_dir'
        _OUTPUT_DIR_STR = 'output'
        _INPUT_DIR_STR = 'input'
        _EXAMPLE_DIR_STR = 'example'

        # output files
        _DONE_FILE = 'done'
        _MEGAHIT_LOG = 'log'

        # full file paths
        self.megahit_path = _MEGAHIT_PATH
        self.unittest_dir = _UNITTEST_DIR_STR
        self.output_dir = os.path.join(_UNITTEST_DIR_STR, _OUTPUT_DIR_STR)
        self.input_dir = os.path.join(_UNITTEST_DIR_STR, _INPUT_DIR_STR)

        # identify the path to megahit directory
        megahit_dir, executable = os.path.split(_MEGAHIT_PATH)

        # find default sample library used for megahit tests
        sample_dir = os.path.join(megahit_dir, _EXAMPLE_DIR_STR)

        # full file paths for each of the interleaved FQ files
        self.sample_1 = os.path.join(sample_dir, _SAMPLE_R1_FQ_STR)
        self.sample_2 = os.path.join(sample_dir, _SAMPLE_R2_FQ_STR)
        self.sample_3 = os.path.join(sample_dir, _SAMPLE_U_FQ_STR)

        # full file paths for each of the trimmed sample files
        self.trim_R1 = os.path.join(self.input_dir, _TRIM_R1_FQ_STR)
        self.trim_R2 = os.path.join(self.input_dir, _TRIM_R2_FQ_STR)
        self.trim_U = os.path.join(self.input_dir, _TRIM_U_FQ_STR)

        # location of generated output files
        self.done_file = os.path.join(self.output_dir, _DONE_FILE)

        # generated output files
        self.megahit_intermed_dir = os.path.join(self.output_dir, _INTERMED_CONTIGS)
        self.megahit_tmp = os.path.join(self.output_dir, _TMP)
        self.megahit_k21 = os.path.join(self.megahit_tmp, _K21)
        self.final_contigs = os.path.join(self.output_dir, _FINAL_CONTIG)
        self.megahit_opts = os.path.join(self.output_dir, _MEGAHIT_OPTS)

        # log files
        self.pipe_log = os.path.join(self.output_dir, _PIPELINE_LOG)
        self.megahit_log = os.path.join(self.output_dir, _MEGAHIT_LOG)

    # TODO: Check if the OSError is thrown in case we remove something improperly
    def clear_dir(self, target_dir):
        """ Selectively remove files in a directory using the given file extension names """

        # output file extensions
        _OUT_EXT = ['.fq', '.fastq', '.fa', '.fasta', '.txt', '.lib', '.bin', '.info', '.lib_info']

        if os.path.exists(target_dir):
            # remove all the files in the intermediate contigs directory
            filelist = [f for f in os.listdir(target_dir) if f.endswith(tuple(_OUT_EXT))]
            for f in filelist:
                f_path = os.path.join(target_dir, f)
                os.remove(f_path)

    # methods to test the formation of the trimmomatic command
    def setUp(self):
        """ create temporary files and directories to be used in the unit tests """

        self.setup_path()
        # create a sample directory to use for input and output
        if not os.path.exists(self.unittest_dir):
            os.makedirs(self.unittest_dir)
        else:
            print("There exists conflicting directory named: {0}".format(self.unittest_dir))

        # create the appropriate directories
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)
        else:
            print("There exists conflicting directory named: {0}".format(self.input_dir))

        if not os.path.exists(self.trim_R1):
            trim_R1_file = open(self.trim_R1, 'w+')
            trim_R1_file.write(_TRIM_R1_FASTQ_STR)
            trim_R1_file.close()
        if not os.path.exists(self.trim_R2):
            trim_R2_file = open(self.trim_R2, 'w+')
            trim_R2_file.write(_TRIM_R2_FASTQ_STR)
            trim_R2_file.close()
        if not os.path.exists(self.trim_U):
            trim_U_file = open(self.trim_U, 'w+')
            trim_U_file.write(_TRIM_U_FASTQ_STR)
            trim_U_file.close()

        input_test_files = [self.sample_1, self.sample_2, self.sample_3]
        for test_file in input_test_files:
            if not os.path.isfile(test_file):
                raise ValueError(
                    "Input file {0} does not exist. Please check megahit/example directory for sample test files".format
                    (test_file))

    def tearDown(self):
        """ delete temporary files and directories generated by setUp method and megahit subprocess calls """

        if os.path.exists(self.unittest_dir):
            if os.path.exists(self.input_dir):
                self.clear_dir(self.input_dir)
                os.rmdir(self.input_dir)
            if os.path.exists(self.output_dir):
                print("the output dir is {0}".format(self.output_dir))
                # first remove all output files
                if os.path.isfile(self.done_file):
                    os.remove(self.done_file)
                if os.path.isfile(self.megahit_log):
                    os.remove(self.megahit_log)
                if os.path.isfile(self.final_contigs):
                    os.remove(self.final_contigs)

                # then remove output subdirectories
                if os.path.exists(self.megahit_intermed_dir):
                    print("the intermediate directory is {0}".format(self.megahit_intermed_dir))
                    self.clear_dir(self.megahit_intermed_dir)
                    os.rmdir(self.megahit_intermed_dir)

                if os.path.exists(self.megahit_tmp):
                    print("the megahit temporary directory is {0}".format(self.megahit_tmp))

                    if os.path.exists(self.megahit_k21):
                        print ("the megahit k21 directory is {0}".format(self.megahit_k21))
                        self.clear_dir(self.megahit_k21)
                        # TODO: WARNING: not sure if this is safe. Would prefer to know what is actually in k21 folder
                        try:
                            os.rmdir(self.megahit_k21)
                        except OSError:
                            print("rmdir called failed. K21 directory removed forcibly with shutil")
                            shutil.rmtree(self.megahit_k21)
                    self.clear_dir(self.megahit_tmp)
                    os.rmdir(self.megahit_tmp)

                self.clear_dir(self.output_dir)
                os.rmdir(self.output_dir)
            # remove the unittest directory
            os.rmdir(self.unittest_dir)
        else:
            print("The unittest directory {0} does not exist".format(self.unittest_dir))

    # region form_trim_cmd_list tests
    def test_form_megahit_cmd_list_no_args(self):
        """ test that form_megahit_cmd_list correctly raises a Value Error when invalid empty string is used in place of
        required in put """

        # arguments to be formatted
        null_megahit_path = ''
        # omit optional arguments
        null_memory = None
        # non-existent fastq files
        null_R1 = ''
        null_R2 = ''
        null_R12 = ''
        null_U = ''
        # option for output directory omitted
        null_outdir = None

        # catch a value exception when strings are empty
        with self.assertRaises(ValueError):
            exec_megahit.form_megahit_cmd_list(null_megahit_path, null_outdir, null_memory, null_R1, null_R2, null_R12,
                                               null_U)

    def test_form_megahit_cmd_list_invalid_num_args(self):
        """ test that form_megahit_cmd_list correctly raises a Type Error when the wrong number of
        input arguments is used """

        # program path
        megahit_path = '/usr/local/megahit/megahit'
        with self.assertRaises(TypeError):
            exec_megahit.form_megahit_cmd_list(megahit_path)

    def test_form_megahit_cmd_list_paired_unpaired_no_max_memory(self):
        """ test that form_megahit_cmd_list correctly generates a megahit command list when passed valid arguments for
        the megahit file path, output directory, paired and unpaired fastq lists """

        # program path
        megahit_path = '/usr/local/megahit/megahit'

        # memory options
        null_memory = None

        # input file names
        test_R1_fq_str = 'left.fastq.gz'
        test_R2_fq_str = 'right.fastq.gz'
        test_U_fq_str = 'unpaired.fastq.gz'

        R1_fq_list = [test_R1_fq_str]
        R2_fq_list = [test_R2_fq_str]
        R12_fq_list = []
        U_fq_list = [test_U_fq_str]

        # no input directory specified
        outdir = '/megahit_unittest_temp_dir/output'

        cmd_no_outdir_list = ['/usr/local/megahit/megahit',
                              '-o', '/megahit_unittest_temp_dir/output',
                              '-1', 'left.fastq.gz',
                              '-2', 'right.fastq.gz',
                              '-r', 'unpaired.fastq.gz']

        self.assertEqual(cmd_no_outdir_list,
                         exec_megahit.form_megahit_cmd_list(megahit_path, outdir, null_memory, R1_fq_list, R2_fq_list,
                                                            R12_fq_list, U_fq_list))

    def test_form_megahit_cmd_list_paired_unpaired(self):
        """ test that form_megahit_cmd_list correctly generates a megahit command list when passed valid arguments for
        the megahit file path, maximum memory, output directory, paired and unpaired fastq files """

        # program path
        megahit_path = '/usr/local/megahit/megahit'

        # memory options
        max_memory = '1'

        # input file names
        test_R1_fq_str = 'left.fastq.gz'
        test_R2_fq_str = 'right.fastq.gz'
        test_U_fq_str = 'unpaired.fastq.gz'

        R1_fq_list = [test_R1_fq_str]
        R2_fq_list = [test_R2_fq_str]
        R12_fq_list = []
        U_fq_list = [test_U_fq_str]

        # no input directory specified
        outdir = '/megahit_unittest_temp_dir/output'

        cmd_no_outdir_list = ['/usr/local/megahit/megahit',
                              '-o', '/megahit_unittest_temp_dir/output',
                              '-m', '1',
                              '-1', 'left.fastq.gz',
                              '-2', 'right.fastq.gz',
                              '-r', 'unpaired.fastq.gz']

        self.assertEqual(cmd_no_outdir_list,
                         exec_megahit.form_megahit_cmd_list(megahit_path, outdir, max_memory, R1_fq_list, R2_fq_list,
                                                            R12_fq_list, U_fq_list))

    def test_form_megahit_cmd_list_interleaved_no_outdir(self):
        """ test that form_megahit_cmd_list correctly generates a megahit commmand list when passed valid arguments for
        the file megahit path, maximum memory, output directory, and interleaved input fastq files """

        # program path
        megahit_path = '/usr/local/megahit/megahit'

        # memory options
        max_memory = '1'

        # input file names
        test_1_fq_str = '/usr/local/megahit/examples/readsInterleaved1.fa.gz'
        test_2_fq_str = '/usr/local/megahit/examples/readsInterleaved2.fa.bz2'
        test_3_fq_str = '/usr/local/megahit/examples/readsInterleaved3.fa'

        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [test_1_fq_str, test_2_fq_str, test_3_fq_str]
        U_fq_list = []

        # no input directory specified
        outdir = '/megahit_unittest_temp_dir/output'

        cmd_no_outdir_list = ['/usr/local/megahit/megahit',
                              '-o', '/megahit_unittest_temp_dir/output',
                              '-m', '1', '--12',
                              '/usr/local/megahit/examples/readsInterleaved1.fa.gz,/usr/local/megahit/examples/readsInterleaved2.fa.bz2,/usr/local/megahit/examples/readsInterleaved3.fa']

        self.assertEqual(cmd_no_outdir_list,
                         exec_megahit.form_megahit_cmd_list(megahit_path, outdir, max_memory, R1_fq_list, R2_fq_list,
                                                            R12_fq_list, U_fq_list))

    def test_form_megahit_cmd_list_paired_unpaired_interleaved(self):
        """ test that form_megahit_cmd_list correctly generates a megahit command list when passed valid arguments for
        the file megahit path, maximum memory, output directory, and paired, unpaired,
        and interleaved input fastq files """

        # program path
        megahit_path = '/usr/local/megahit/megahit'

        # memory options
        max_memory = '1'

        # input file names
        test_1_fq_str = '/usr/local/megahit/examples/readsInterleaved1.fa.gz'
        test_2_fq_str = '/usr/local/megahit/examples/readsInterleaved2.fa.bz2'
        test_3_fq_str = '/usr/local/megahit/examples/readsInterleaved3.fa'

        # input file names
        test_R1_fq_str = '/usr/local/megahit/examples/left.fastq.gz'
        test_R2_fq_str = '/usr/local/megahit/examples/right.fastq.gz'
        test_U_fq_str = '/usr/local/megahit/examples/unpaired.fastq.gz'

        R12_fq_list = [test_1_fq_str, test_2_fq_str, test_3_fq_str]
        R1_fq_list = [test_R1_fq_str]
        R2_fq_list = [test_R2_fq_str]
        U_fq_list = [test_U_fq_str]

        # no input directory specified
        outdir = '/megahit_unittest_temp_dir/output'

        # full command to match
        cmd_no_outdir_list = ['/usr/local/megahit/megahit',
                              '-o', '/megahit_unittest_temp_dir/output',
                              '-m', '1',
                              '-1', '/usr/local/megahit/examples/left.fastq.gz',
                              '-2', '/usr/local/megahit/examples/right.fastq.gz',
                              '--12',
                              '/usr/local/megahit/examples/readsInterleaved1.fa.gz,/usr/local/megahit/examples/readsInterleaved2.fa.bz2,/usr/local/megahit/examples/readsInterleaved3.fa',
                              '-r', '/usr/local/megahit/examples/unpaired.fastq.gz']

        self.assertEqual(cmd_no_outdir_list,
                         exec_megahit.form_megahit_cmd_list(megahit_path, outdir, max_memory, R1_fq_list, R2_fq_list,
                                                            R12_fq_list, U_fq_list))
    # endregion

    # region run_megahit_cmd_list tests
    def test_run_megahit_cmd_list_no_args(self):
        """ test shall check that run_megahit correctly raises a Value Error when invalid empty string is used in place
        of required input """

        # arguments to be formatted
        null_megahit_path = ''
        # omit optional arguments
        null_memory = None
        # non-existant fastq files
        null_R1 = ''
        null_R2 = ''
        null_R12 = ''
        null_U = ''
        # option to omit output directory
        null_outdir = None
        # option to omit output log
        null_log = None

        with self.assertRaises(ValueError):
            exec_megahit.run_megahit(null_megahit_path, null_outdir, null_memory, null_R1, null_R2, null_R12, null_U)

    def test_run_megahit_cmd_list_invalid_num_args(self):
        """ test that run_megahit correctly raises a Type Error when an invalid number of arguments is used in place of
         required input """

        # program path
        megahit_path = '/usr/local/megahit/megahit'
        with self.assertRaises(TypeError):
            exec_megahit.run_megahit(megahit_path)

    def test_run_megahit_cmd_list_interleaved_good_stderr(self):
        """ test shall check if megahit subproccess call reports successful completion message to stderr when
        run_megahit is passed valid arguments for the megahit file path, output directory, maximum memory, and
        interleaved fastq files """

        # indicators of successful megahit execution
        _SUCCESSFUL_COMPLETION = 'ALL DONE'

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [self.sample_1, self.sample_2, self.sample_3]
        U_fq_list = []

        # no input directory specified
        null_log = None

        _, stderr = exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list,
                                             R12_fq_list, U_fq_list)

        self.assertNotEqual(stderr.find(_SUCCESSFUL_COMPLETION), -1)


    def test_run_megahit_bad_input_cmd_list(self):
        """ test that megahit subproccess call does not report successful completion message to stderr when
        run_megahit is not passed valid passing arguments for the paired and unpaired input fastq files """

        # indicators of successful megahit execution
        _SUCCESSFUL_COMPLETION = 'ALL DONE'

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = [self.trim_R1]
        R2_fq_list = [self.trim_R2]
        R12_fq_list = []
        U_fq_list = [self.trim_U]

        # no input directory specified
        null_log = None

        _, stderr = exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list,
                                             R12_fq_list, U_fq_list)

        self.assertEqual(stderr.find(_SUCCESSFUL_COMPLETION), -1)

    def test_run_megahit_cmd_list_interleaved_assembly_occurred(self):
        """ test that megahit log is generated when run_megahit is passed valid arguments for the megahit file path,
        output directory, maximum memory, and interleaved fastq input files, indicating modified input """

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [self.sample_1, self.sample_2, self.sample_3]
        U_fq_list = []

        # no input directory specified
        null_log = None

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list, R12_fq_list,
                                 U_fq_list)
        # check that an assembly was started
        self.assertTrue(os.stat(self.megahit_log).st_size > 0)

    def test_run_megahit_cmd_list_interleaved_completion_file_generated(self):
        """ test that done file is generated when run_megahit is passed valid arguments for the megahit file path,
        output directory, maximum memory, and interleaved fastq files, indicating successful completion """

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [self.sample_1, self.sample_2, self.sample_3]
        U_fq_list = []

        # no input directory specified
        null_log = None

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list, R12_fq_list,
                                 U_fq_list)
        # check that an assembly was started
        self.assertTrue(os.path.exists(self.done_file))

    def test_run_megahit_cmd_list_trimmed_completion_file_not_generated(self):
        """ test that done file is not generated when run_megahit is passed invalid input fastq files """

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = [self.trim_R1]
        R2_fq_list = [self.trim_R2]
        R12_fq_list = []
        U_fq_list = [self.trim_U]

        # no input directory specified
        null_log = None

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list, R12_fq_list,
                                 U_fq_list)

        self.assertFalse(os.path.exists(self.done_file))

    def test_run_megahit_cmd_list_interleaved_final_contigs_generated(self):
        """ test that a final contig.fasta file is generated when run_megahit is passed valid arguments for the
        megahit file path, output directory, maximum memory, and interleaved fastq files,
        indicating successful completion """

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [self.sample_1, self.sample_2, self.sample_3]
        U_fq_list = []

        # no input directory specified
        null_log = None

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list, R12_fq_list,
                                 U_fq_list)
        # check that an assembly was started
        self.assertTrue(os.path.exists(self.final_contigs))

    def test_run_megahit_cmd_list_trimmed_final_contigs_not_generated(self):
        """ test that final contig.fasta file is not generated when run_megahit is passed invalid arguments
         for the input fastq files """

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = [self.trim_R1]
        R2_fq_list = [self.trim_R2]
        R12_fq_list = []
        U_fq_list = [self.trim_U]

        # no input directory specified
        null_log = None

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list, R12_fq_list,
                                 U_fq_list)

        self.assertFalse(os.path.exists(self.final_contigs))

    # TODO: future tests to incorporate
    """
    def test_run_megahit_cmd_list_interleaved_log_generated(self):
        logging.basicConfig(filename=self.pipe_log, level=logging.INFO)
        logging.debug('unittest: test_run_megahit_cmd_list_interleaved_log_generated')

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = []
        R2_fq_list = []
        R12_fq_list = [self.sample_1, self.sample_2, self.sample_3]
        U_fq_list = []

        exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list, R2_fq_list,
                                 R12_fq_list, U_fq_list, self.pipe_log)
        # check that an assembly was started
        self.assertTrue(os.stat(self.pipe_log).st_size > 0)
    """
    '''
    def test_run_megahit_cmd_list_fq_format(self):
        # indicators of successful megahit execution
        _SUCCESSFUL_COMPLETION = 'ALL DONE'

        # memory options
        max_memory = '.8'

        # list of fq and fa files to be passed to megahit
        R1_fq_list = [self.sample_1]
        R2_fq_list = [self.sample_2]
        R12_fq_list = []
        U_fq_list = [self.sample_3]

        # no input directory specified
        null_log = None

        stdout, stderr = exec_megahit.run_megahit(self.megahit_path, self.output_dir, max_memory, R1_fq_list,
                                                  R2_fq_list,
                                                  R12_fq_list, U_fq_list, null_log)

        print(stderr)

        self.assertEqual(stderr.find(_SUCCESSFUL_COMPLETION), -1)
    '''

