""" Module to test exec_trim member methods. This module contains unit tests for exec_trim.py. """

import os
import unittest
# import logging

from src import exec_trim

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"

# input files contents
_SAMPLE_R1_FASTQ_STR = """@K00180:79:H3NG2BBXX:2:1101:28199:1086 1:N:0:NCCAAT+NGATCTCG
NGAAAAATTTCCTATAGCTTACTTGCAAGGTGCATGATTGAAACGGCAACTATTTATGAGACTTAAGAATTTCTTGACAGAACTGCAGCTAAAATATAGCTTGCGGGTACGGAAAACTGCGGAGCATGCTTCAATAGTGGTATCAGCTTC
+
#AAFFJJFJJJJJJAJJJ7JJJJJJJJJJJJJJJ-FJJJFJJJJJJJJJJAJJJJJAJJ<JJJJJJJJJJJJJJFJFFJFJJJJJJJ<JJJJJJJJ<JJJJJFF<FJAJF<FJJJJJ-7AA)FAJJFFJJFJ<<-FFAJ--<AAA7<<FF
@K00180:79:H3NG2BBXX:2:1101:28666:1086 1:N:0:NCCAAT+NGATCTCG
NAACGTTATCAGCAAGCAATACCCGCTACAGATCGATATCTATACCACCCATGAGGACAAACTAATCATATCCAATAAAATACATCCTAAGTCGGAAGATAGCAAAGATAGCGGTATCTGATTGAAGAACCTATGGGGGCGCTACCGGAT
+
#AAFF<FJFJJFJJJJJJJJJJJJJJJJJJJJJJJJFJFJ<JFJJJJJ<JFFFJFAJJJJJJFJJJAFAJ-77JJFFJJJFJJJAFJFAJ<AJJAJJJJJJFFJJJJ-7AFFAF)7<A-7<7--AJJJJJJ7F<<A<A--)--AF)AAFA
@K00180:79:H3NG2BBXX:2:1101:28828:1086 1:N:0:NCCAAT+NGATCTCG
NTGAAGACTACCTGACTTTGAGGTGAGAAGAATGAAAGCTAAAAACCATAATACAAACAGAGAAAAACGTGTCATGAAAAAACGGTTCAGGACAGCAGTCTTTCTTCTTCTTTTTGCGGGAAGTGGTTTGACTGTTTTTCGATATTTCAA
+
#A<FFJJJJJJJJJFJJJ7AJJJAJJJJJJFJJJJJJJFFFFJJJJFF-F<F<JJJJJJJJAFJJJJJJFJ-AFFAFFJJJJJFJ7AF-AAJJFJJJFAAAFJJJJJJFA7FFJA7F-77)-)7-A7-A-<A-<<<F<FA-7-7-<AFF-
@K00180:79:H3NG2BBXX:2:1101:28970:1086 1:N:0:NCCAAT+NGATCTCG
NAGGTCGAGCAGGATGCCTGCCAGCACCTTCTTGTAACGGATCGTATCGCCGCGCAGCGCGTCGGTAGTGACGTTGTCGATGTTCTTGATGAACACCACGTCGGCGTCGATCTCGGTCAGGTGTTCGATCAGCGCGCCGTGGCCCGCAGG
+
#AAA<FJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJ<JF<JA<J<J-FFAJJJJ7FAJJJJJ<AJJ-FAJFJ7FAA<AFFAFFAJFFJ<FJJFJJJAAJF77F<-A-F-A-AF)A<AJ<-)7<<AA7A<<A<<7<F)7-7<<F-A77F
@K00180:79:H3NG2BBXX:2:1101:29112:1086 1:N:0:NCCAAT+NGATCTCG
NATTTATTCAGAATATTTAAATAAAATATTCACTTCCAAGTATTGGAAAGATGAGATCCACAAATATAGCATGGGAACAACCGTTCACACATTGACAATAAGTAGAGCAAATGAGGTTGTAATCCCTTTTGCTAGCGAACAAGATCAAAT
+
#AAAFFJJJJFJJJJFFJJJJAJJJJJJJFJJJF-F-7<<<F<JA7FFJ<FFJJ<<FFJJJJJJFJ-A<FJ-77AFAFJJJJ7F-AFJJJJ<FFFFFJ7FJAAF<FF-F-F--AF7---7A7A-777--<A<FF-))7<AFJ<<--<<<-"""

_SAMPLE_R2_FASTQ_STR = """@K00180:79:H3NG2BBXX:2:1101:28199:1086 2:N:0:NCCAAT+NGATCTCG
NACCTACCCTTAACGTCTCTGCCAAGTGCTGATAAGAAGAAGCATCCCCGCCACCAACTATATGTAACTTCACTGAACAGTTATTCATTTGTGCCACTGTGCGGATAGCTAATGTCAGCTGTATACGAAAATTCAACTTTCCTACCCAAC
+
#A<<FFJ7FF7FJFF-F<J77FAJA-<AF7-AAAF-7-A<AFJJAJJFJJFJJJJJFFFJFJJ<FAJJFFJFJ<JFA-7F--77-<<--<-A-7)7A-7A-A)--7AFA-7<<--77-<A)7--77-A---7-7<A<-77A<7--)7--)
@K00180:79:H3NG2BBXX:2:1101:28666:1086 2:N:0:NCCAAT+NGATCTCG
NCACCTGCTCGAATATCTTGAACGCCGATCCGTCCGCTAAATGGATATCCATGAAGATGAGATCCGGATGAGGCGATTTCCTGAAGAAATCCACGCTCTCCTCGATACTCTCCAACTCGGCGAGCACTTCGAACTCAGTAACGCAACTCT
+
#7<<F-FJAFJFFFJFJ<JJJAJJJFAA-7FAJJ7FAAFF<7FJJJ<-FFJA<-FAJ-7F7FJF7F-<A7<-<FFA7--AJ-<--AJ-F-7F<F)A7FF<A-<7)-7<-77<)7-A-<F)-7)-7<-F-<--77--7<7------7)---
@K00180:79:H3NG2BBXX:2:1101:28828:1086 2:N:0:NCCAAT+NGATCTCG
NCAGTCGTTACGATATCATTCCCGTTTTGTATCTCATCTTCTAAATTTTCCTGTAAATCAATATATGCCGTTTCTCCGGTCAGCATTTTGTAGTTACCGTGCGTTGACAGAAAATAAAAATTTAAAAATCAGGAGACGTTTTGAATGTAC
+
#A<-<7FA-AJ<FJ7-<-<7F-7F7A-AJA7-AJJFAJJJ7-<-A-AFAJFF-FJJ--7F7-7-7A<A77----A<-)))-A))-7-<--J---7<<-<<)))<-)-7<--7---77<7<----7--7A-----)-))-----------7
@K00180:79:H3NG2BBXX:2:1101:28970:1086 2:N:0:NCCAAT+NGATCTCG
NCCGTCACCTCCCACGCCTATCCCGAAGGCGCGCGCAAGGCCGTCGAGGAACATCTGGTCGAAGGGGCGGTCTATGCGGCCGCGAACGGCGTGGCGAAGATTCATTTCACGGTTGCTCCCGAACACATCGAAGGCGTTCAGAAGCGTCTC
+
#-<<<AA<J<-7JJJAJA7F<AFJ7A<777-77J-7-<-FF-7-7-777-AFA7AA--7AFF77--7A-77A-A-77-7AJAJA<AJF)<--)-<-7<7<--<FA-AFJF-))-)7-7<J)77AJ-77A----7)-))7---7--7-7)7
@K00180:79:H3NG2BBXX:2:1101:29112:1086 2:N:0:NCCAAT+NGATCTCG
NCATTACCTCTAATTGCAGGAGGGTTGGCGTAAGCGCATCGATACTGGCTTCTGTGGTAAACAGCGTGTCGAAATGCTGACGGATTCGCACCCAGCTGTCACAAAGCTCTCTGGTATTTTTGCTGTCCTTACGCGTGGCACGGCGGGTCT
+
#7-<<FF7--<---<F--<7A-7----77--7FF7<-7AA----<7-AF--7-7-777AJJFJ<A7---7--7A--7--AF)A<A<77-7))<)-))77<<--7<AA<-7)-7)7--7A<-<7---7)-7----)))-<-<<))<)-)-7"""

# message printed when trimmomatic process completes successfully
_SUCCESSFUL_COMPLETION = 'Completed successfully'

# TODO: CHANGE CONSTANT ACCORDING TO SYSTEM INSTALL
_TRIMMOMATIC_PATH = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'


class TestExecTrim(unittest.TestCase):
    """This module shall act as a framework for exec_trim unit tests """

    @staticmethod
    def _print_msg_re_existing(what_exists, exists_name):
        """ static method for reporting existing files and directories set up for each unit test

        Args:
            what_exists: the type of file which is reported to exist
            exists_name: the name of the file
        """

        print("A {0} named '{1}' already exists, and will be used.".format(what_exists, exists_name))

    def _setup_path(self):
        """ creates strings corresponding to the temporary directories and files to be used in the unit tests """

        # file and directory names to be used as output files
        _FQ_R1_STR = 'sample_R1_20.fastq'
        _FQ_R2_STR = 'sample_R2_20.fastq'

        # log files
        _PIPELINE_LOG = 'pipeline.log'
        _TRIM_LOG = 'trim.log'

        # temporary directories to be used
        _UNITTEST_DIR_STR = 'trimmomatic_unittest_temp_dir'
        _INPUT_DIR_STR = 'input'
        _OUTPUT_DIR_STR = 'output'

        # full file paths
        self.trim_path = _TRIMMOMATIC_PATH
        self.unittest_dir = _UNITTEST_DIR_STR
        self.input_dir = os.path.join(_UNITTEST_DIR_STR, _INPUT_DIR_STR)
        self.output_dir = os.path.join(_UNITTEST_DIR_STR, _OUTPUT_DIR_STR)
        self.fq_R1 = os.path.join(self.input_dir, _FQ_R1_STR)
        self.fq_R2 = os.path.join(self.input_dir, _FQ_R2_STR)
        self.trim_log = os.path.join(self.output_dir, _TRIM_LOG)
        self.pipe_log = os.path.join(self.unittest_dir, _PIPELINE_LOG)

    # region form_trim_cmd_list tests
    def setUp(self):
        """ creates temporary directories and files to be used in the unit tests """

        self._setup_path()
        # create temporary test directories
        temp_dir_list = [self.unittest_dir, self.input_dir, self.output_dir]
        for temp_dir in temp_dir_list:
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            else:
                self._print_msg_re_existing("directory", temp_dir)
                # create a sample input paired end fastq files
        if not os.path.isfile(self.fq_R1):
            fq_R1_file = open(self.fq_R1, "w+")
            fq_R1_file.write(_SAMPLE_R1_FASTQ_STR)
            fq_R1_file.close()
        else:
            self._print_msg_re_existing("fastq file", self.fq_R1)
        if not os.path.isfile(self.fq_R2):
            fq_R2_file = open(self.fq_R2, "w+")
            fq_R2_file.write(_SAMPLE_R2_FASTQ_STR)
            fq_R2_file.close()
        else:
            self._print_msg_re_existing("fastq file", self.fq_R2)

    def tearDown(self):
        """ deletes temporary directories and files generated by setup and also by the trimmomatic proccess call """

        # output file extensions
        _OUT_EXT = [".fq", ".fastq", ".log"]

        if os.path.exists(self.unittest_dir):
            # first remove all output files
            if os.path.isfile(self.fq_R1):
                os.remove(self.fq_R1)
            if os.path.isfile(self.fq_R2):
                os.remove(self.fq_R2)
            if os.path.exists(self.pipe_log):
                os.remove(self.pipe_log)

            # then remove input and output subdirectories
            dirs_to_remove = [self.input_dir, self.output_dir]
            for curr_dir in dirs_to_remove:
                if os.path.exists(curr_dir):
                    file_list = [f for f in os.listdir(curr_dir) if f.endswith(tuple(_OUT_EXT))]
                    for f in file_list:
                        f_path = os.path.join(curr_dir, f)
                        os.remove(f_path)
                    os.rmdir(curr_dir)
            # remove the unittest directory
            os.rmdir(self.unittest_dir)
        else:
            print("The unittest directory {0} does not exist".format(self.unittest_dir))

    # region for the testing of form_trim_cmd_list
    def test_form_trim_cmd_list_no_args(self):
        """ test that form_trim_cmd_list raises a Value Error when invalid empty string is used in place of
        required input """

        # arguments to be formatted
        null_trimmomatic_path = ''
        no_mode = ''
        null_fq = ''
        null_adapter_opt = None
        # option for output directory omitted
        null_outdir = None
        null_trim_log = None

        with self.assertRaises(ValueError):
            output, err = exec_trim.form_trim_cmd_list(null_trimmomatic_path, no_mode, null_fq, null_adapter_opt,
                                                       null_outdir, null_trim_log)

    def test_form_trim_cmd_list_invalid_mode(self):
        """ test that from_trim_cmd_list raises a Value Error when user fails to specify paired or single ended mode"""

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'invalid'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        null_outdir = None
        # option to retain trim log not specified
        null_trim_log = None

        with self.assertRaises(ValueError):
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, null_outdir, null_trim_log)

    def test_form_trim_cmd_list_no_adapter_opt(self):
        """ test that form_trim_cmd_list raises a Value error when user fails to specify trimmomatic clip settings """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'PE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        null_adapter_opt = []
        # option to redirect output not specified
        null_outdir = None
        # option to retain trim log not specified
        null_trim_log = None

        with self.assertRaises(ValueError):
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, null_adapter_opt, null_outdir,
                                         null_trim_log)

    def test_form_trim_cmd_list_SE_wrong_number_fastq_input(self):
        """ test that form_trim_cmd_list rasies a Value Error if more than a single file is specified for
        single-end mode """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'SE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        null_outdir = None
        # option to retain trim log not specified
        null_trim_log = None

        with self.assertRaises(ValueError):
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, null_outdir, null_trim_log)

    def test_form_trim_cmd_list_PE_wrong_number_fastq_input(self):
        """ test that form_trim_cmd_list correctly raises a Value Error if more than two fastq files are specified for
        paired-end mode """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'PE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        fq2_R1 = 'trimmomatic_unittest_temp_dir/input/sample2_R1_20.fastq'
        fq2_R2 = 'trimmomatic_unittest_temp_dir/input/sample2_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2, fq2_R1, fq2_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        null_outdir = None
        # option to retain trim log not specified
        null_trim_log = None

        with self.assertRaises(ValueError):
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, null_outdir, null_trim_log)

    def test_form_trim_cmd_list_PE_no_outdir_no_trim_log(self):
        """ test that form_trim_cmd_list correctly generates a paired-end trimmomatic command without a trim log when
        passed valid arguments for the trimmomatic file path, input fastq """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'PE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        null_outdir = None
        # option to retain trim log not specified
        null_trim_log = None

        cmd_no_outdir_no_trim_log_list = ['java', '-jar', '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar', 'PE',
                                          'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq',
                                          'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq',
                                          '-baseout', 'trimmomatic_unittest_temp_dir/input/sample_20.fastq',
                                          'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:15', 'MINLEN:70']

        self.assertEqual(exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, null_outdir,
                                                      null_trim_log), cmd_no_outdir_no_trim_log_list)

    def test_form_trim_cmd_list_PE_no_trim_log(self):
        """ test that form_trim_cmd_list correctly generates a paired-end trimmomatic command without a trim log
        when passed valid arguments for the trimmomatic file path, input fastq, and output directory """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'PE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        outdir = 'trimmomatic_unittest_temp_dir/output/'
        # option to retain trim log not specified
        null_trim_log = None

        cmd_no_outdir_no_trim_log_list = ['java', '-jar', '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar', 'PE',
                                          'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq',
                                          'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq',
                                          '-baseout', 'trimmomatic_unittest_temp_dir/output/sample_20.fastq',
                                          'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:15', 'MINLEN:70']

        self.assertEqual(
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, outdir, null_trim_log), cmd_no_outdir_no_trim_log_list)

    def test_form_trim_cmd_list_PE(self):
        """ test that form_trim_cmd_list correctly generates a paired-end trimmomatic command when passed valid
        arguments for the trimmomatic file path, inmput fastq, and ouput directory """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'PE'
        # paired fq input
        fq_R1 = 'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq'
        fq_R2 = 'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq'
        input_fq_list = [fq_R1, fq_R2]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        outdir = 'trimmomatic_unittest_temp_dir/output/'
        # option to retain trim log not specified
        trim_log = 'trimmomatic_unittest_temp_dir/output/trim.log'

        cmd_no_outdir_no_trim_log_list = ['java', '-jar', '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar', 'PE',
                                          'trimmomatic_unittest_temp_dir/input/sample_R1_20.fastq',
                                          'trimmomatic_unittest_temp_dir/input/sample_R2_20.fastq',
                                          '-baseout', 'trimmomatic_unittest_temp_dir/output/sample_20.fastq',
                                          '-trimlog', 'trimmomatic_unittest_temp_dir/output/trim.log',
                                          'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:15', 'MINLEN:70']

        self.assertEqual(
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, outdir, trim_log),
            cmd_no_outdir_no_trim_log_list)

    def test_form_trim_cmd_list_SE(self):
        """ test that form_trim_cmd_list correctly generates a single-end trimmomatic command when passed valid
        arguments for the trimmomatic file path, input fastq, and output directory """

        # program path
        trimmomatic_path = '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar'
        # paired end trimming
        mode = 'SE'
        # paired fq input
        fq_unpaired = 'trimmomatic_unittest_temp_dir/input/sample_20.fastq'
        input_fq_list = [fq_unpaired]
        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")
        # option to redirect output not specified
        outdir = 'trimmomatic_unittest_temp_dir/output/'
        # option to retain trim log not specified
        trim_log = 'trimmomatic_unittest_temp_dir/output/trim.log'

        cmd_no_outdir_no_trim_log_list = ['java', '-jar', '/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar', 'SE',
                                          'trimmomatic_unittest_temp_dir/input/sample_20.fastq',
                                          'trimmomatic_unittest_temp_dir/output/sample_20.fastq',
                                          '-trimlog', 'trimmomatic_unittest_temp_dir/output/trim.log',
                                          'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:15', 'MINLEN:70']

        self.assertEqual(
            exec_trim.form_trim_cmd_list(trimmomatic_path, mode, input_fq_list, adapter_opt, outdir, trim_log),
            cmd_no_outdir_no_trim_log_list)
    # endregion

    # region to test run_trimmomatic
    def test_run_trimmomatic_PE_good_stderr(self):
        """ test that trimmomatic subprocess call reports successful completion to stderr when run_trimmomatic is passed
        valid arguments for the trimmomatic file path, trim mode, input fastq files, adapter options, and trim log """

        # input fastq paths
        input_fq_list = [self.fq_R1, self.fq_R2]

        # paired end trimming
        mode = 'PE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        output, err = exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir,
                                                self.trim_log)

        self.assertNotEqual(err.find(_SUCCESSFUL_COMPLETION), -1)

    def test_run_trimmomatic_PE_output_exists(self):
        """ test that a set of trimmed fastq files are generated in the form of paired forward, paired reverse,
        unpaired forward, and unpaired reverse output files when run_trimmomatic is passed valid arguments for
        the trimmomatic file path, trim mode, input fastq files, adapter options, and trim log """

        fq_out_1U_str = 'sample_20_1U.fastq'
        fq_out_1P_str = 'sample_20_1P.fastq'
        fq_out_2U_str = 'sample_20_2U.fastq'
        fq_out_2P_str = 'sample_20_2P.fastq'

        # input fastq paths
        input_fq_list = [self.fq_R1, self.fq_R2]
        fq_out_1U_str = os.path.join(self.output_dir, fq_out_1U_str)
        fq_out_1P_str = os.path.join(self.output_dir, fq_out_1P_str)
        fq_out_2U_str = os.path.join(self.output_dir, fq_out_2U_str)
        fq_out_2P_str = os.path.join(self.output_dir, fq_out_2P_str)

        # paired end trimming
        mode = 'PE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir, self.trim_log)

        output_cnt = 0

        if os.path.isfile(fq_out_1U_str) > 0:
            output_cnt += 1
        if os.path.isfile(fq_out_1P_str) > 0:
            output_cnt += 1
        if os.path.isfile(fq_out_2U_str) > 0:
            output_cnt += 1
        if os.path.isfile(fq_out_2P_str) > 0:
            output_cnt += 1

        self.assertEqual(output_cnt, 4)

    def test_run_trimmomatic_PE_fastq_modified(self):
        """ test that trim log is generated when run_trimmomatic is passed valid arguments for the trimmomatic file path,
        trim mode, input fastq files, adapter options, and trim log, indicating non-zero changes made by trimmomatic on
        the paired end input fastq files """

        # input fastq paths
        input_fq_list = [self.fq_R1, self.fq_R2]

        # paired end trimming
        mode = 'PE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir, self.trim_log)
        self.assertTrue(os.stat(self.trim_log).st_size > 0)

    def test_run_trimmomatic_SE_output_exists(self):
        """ test that a trimmed fastq file is generated when run_trimmomatic is passed valid arguments for the
        trimmomatic file path, trim mode, input fastq files, adapter options, and trim log """

        # temporary directories to be used
        fq_out_str = 'sample_20.fastq'

        # input fastq paths
        input_fq_list = [self.fq_R1]

        # output fastq path
        fq_out_str = os.path.join(self.output_dir, fq_out_str)

        # paired end trimming
        mode = 'SE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir, self.trim_log)
        self.assertTrue(os.path.isfile(fq_out_str))

    def test_run_trimmomatic_SE_output_modified(self):
        """ test that trim log is generated when run_trimmomatic is passed valid arguments for the trimmomatic file path,
        trim mode, input fastq files, adapter options, and trim log, indicating non-zero changes made by trimmomatic
        on the single end input fastq files """

        # input fastq paths
        input_fq_list = [self.fq_R1]

        # paired end trimming
        mode = 'SE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir, self.trim_log)

        self.assertTrue(os.stat(self.trim_log).st_size > 0)

    '''
    TODO: Not concerned with unit testing the log writing capabilites at the moment
    TODO: Test passes when run individually but passes when run with other commands
    def test_run_trimmomatic_PE_pipe_log_exists(self):
        logging.basicConfig(filename=self.pipe_log, level=logging.INFO)
        logging.debug("unittest: test_run_trimmomatic_PE_pipe_log_exists")

        # input fastq paths
        input_fq_list = [self.fq_R1, self.fq_R2]

        # paired end trimming
        mode = 'PE'

        # adapter options
        adapter_opt_line = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70"
        adapter_opt = adapter_opt_line.split(" ")

        exec_trim.run_trimmomatic(self.trim_path, mode, input_fq_list, adapter_opt, self.output_dir,
                                  self.trim_log, self.pipe_log)

        print( self.pipe_log )
        self.assertTrue(os.stat(self.pipe_log).st_size > 0)

    '''
