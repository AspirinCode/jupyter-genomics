""" Module to test exec_fastqc member methods. This module contains unit tests for exec_fatqc.py"""

import os
import unittest

from src import exec_fastqc

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"

# input file contents
_SHORT_FASTQ_STR = """@K00180:79:H3NG2BBXX:1:1115:6888:46135 1:N:0:GCCAAT+AGATCTCG
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
NAGGTCGAGCAGGATGCCTGCCAGCACCTTCTTGTAACGGATCGTATCGCCGCGCAGCGCGTCGGTAGTGACGTTGTCGATGTTCTTGATGAACACCACGTCGGCGTCGATCTCGGTCAGGTGTTCGATCAGCGCGCCGTGGCCCGCAGG
+
#AAA<FJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJ<JF<JA<J<J-FFAJJJJ7FAJJJJJ<AJJ-FAJFJ7FAA<AFFAFFAJFFJ<FJJFJJJAAJF77F<-A-F-A-AF)A<AJ<-)7<<AA7A<<A<<7<F)7-7<<F-A77F
@K00180:79:H3NG2BBXX:2:2225:10247:10950 1:N:0:GCCAAT+AGATCTCG
GCCATAGCCGGGATATCCTGATACCGCATATCGCCGACTGGCAGCCGTCTGGCTTCCATCCGCTCATCCCAAAAAGGCGTCCCTTCTTGGCTGCGGCAGCATACGTAGAACCTATGTCGTCGTCGGAACAAAGTCCGTTCCACTCCGCCC
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJ<JFA-7AJJJJJFJJFJJJJA-AJ7JJJJFAFFJJJJJJFFJJJFJJJJJAFFJFJJJ7F<-AA7<)F<JAFJJJAFAJ77-F-A7AA7F--"""

# TODO: CHANGE DEPENDING ON SYSTEM INSTALLATION
_FASTQC_PATH = '/usr/local/FastQC/./fastqc'


class TestExecFastQC(unittest.TestCase):
    """Unit test exec_fastqc methods """

    # output file extensions
    _OUT_EXT = [".html", ".zip", ".log"]

    @staticmethod
    def _print_msg_re_existing(what_exists, exists_name):
        """Print message about already-existing files and directories.

        Args:
            what_exists(str): the type of object (file, directory, etc) that has been found to already exist
            exists_name(str): the name of the object that has been found to already exist
        """

        print("A {0} named '{1}' already exists, and will be used.".format(what_exists, exists_name))

    def _setup_paths(self):
        """Create strings corresponding to the temporary directories and files to be used in the unit tests."""

        # file and directory names to be used as output files
        _SHORT_FQ_STR = 'short_sample.fastq'
        _LONG_FQ_STR = 'long_sample.fastq'

        # log file
        _LOG_STR = 'out.log'

        # temporary directories to be used
        _UNITTEST_DIR_STR = 'fastqc_unittest_temp_dir'
        _INPUT_DIR_STR = 'input'
        _OUTPUT_DIR_STR = 'output'

        self.fastqc_path = _FASTQC_PATH
        self.unittest_dir = _UNITTEST_DIR_STR
        self.input_dir = os.path.join(_UNITTEST_DIR_STR, _INPUT_DIR_STR)
        self.output_dir = os.path.join(_UNITTEST_DIR_STR, _OUTPUT_DIR_STR)
        self.short_fq = os.path.join(self.input_dir, _SHORT_FQ_STR)
        self.long_fq = os.path.join(self.input_dir, _LONG_FQ_STR)
        self.log_file = os.path.join(self.unittest_dir, _LOG_STR)

    # TODO: Check if the OSError is thrown in case we remove something improperly
    def clear_dir(self, target_dir):
        """Remove all html, zip, and/or log files in the given directory."""

        if os.path.exists(target_dir):
            # remove all the files in the intermediate contigs directory
            filelist = [f for f in os.listdir(target_dir) if f.endswith(tuple(self._OUT_EXT))]
            for f in filelist:
                f_path = os.path.join(target_dir, f)
                os.remove(f_path)

    def setUp(self):
        """ create temporary directories and files to be used in the unit tests """

        # used to create the longer fastq file
        _FASTQC_CYCLE = 50

        # create path names
        self._setup_paths()

        # create directories for input and output
        temp_dir_list = [self.unittest_dir, self.input_dir, self.output_dir]
        for temp_dir in temp_dir_list:
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
            # do not overwrite previous files
            else:
                self._print_msg_re_existing("directory", temp_dir)

        # generate input fastq files for various test cases
        if not os.path.isfile(self.short_fq):
            short_fq_file = open(self.short_fq, "w+")
            short_fq_file.write(_SHORT_FASTQ_STR)
            short_fq_file.close()
        else:
            self._print_msg_re_existing("fastq file", self.short_fq)

        if not os.path.isfile(self.long_fq):
            long_fq_file = open(self.long_fq, "a")
            for count in range(0, _FASTQC_CYCLE):
                long_fq_file.write(_SHORT_FASTQ_STR)
            long_fq_file.close()
        else:
            self._print_msg_re_existing("fastq file", self.long_fq)

    def tearDown(self):
        """ delete temporary directories and files generated by setup and also by the fastqc proccess call """

        if os.path.exists(self.unittest_dir):
            # first remove all output files
            if os.path.isfile(self.short_fq):
                os.remove(self.short_fq)
            if os.path.isfile(self.long_fq):
                os.remove(self.long_fq)

            # then remove input and output subdirectories
            if os.path.isfile(self.log_file):
                print("removing: {0}".format(self.log_file))
                os.remove(self.log_file)
            self.clear_dir(self.unittest_dir)

            dirs_to_remove = [self.input_dir, self.output_dir]
            for curr_dir in dirs_to_remove:
                if os.path.exists(curr_dir):
                    file_list = [f for f in os.listdir(curr_dir) if f.endswith(tuple(self._OUT_EXT))]
                    for f in file_list:
                        f_path = os.path.join(curr_dir, f)
                        os.remove(f_path)
                    os.rmdir(curr_dir)
            # remove the unittest directory
            os.rmdir(self.unittest_dir)
        else:
            print("The temporary unittest directory {0} does not exist".format(self.unittest_dir))

    # region form_fastqc_cmd_list tests
    def test_form_fastqc_cmd_list_no_args(self):
        """ test that form_fastqc_cmd_list correctly raises a Value Error when invalid empty string is used in place
        of required input """

        # arguments to be formatted
        null_fq_path = ''
        null_fq = ''
        # option for output directory omitted
        null_outdir = None

        # catch a value exception when strings are empty
        with self.assertRaises(ValueError):
            exec_fastqc.form_fastqc_cmd_list(null_fq_path, null_fq, null_outdir)

    def test_form_fastqc_cmd_list_no_outdir(self):
        """ test that form_fastq_cmd_list correctly generates a fastqc command list when passed valid arguments for
        the fastqc path and fastq input files """

        # arguments to be formatted
        fastqc_path = '/usr/local/FastQC/./fastqc'
        fq = 'sample.fq'
        # option to redirect output not specified
        null_outdir = None

        cmd_no_outdir_list = ['/usr/local/FastQC/./fastqc', 'sample.fq']

        # noinspection PyTypeChecker
        self.assertEqual(exec_fastqc.form_fastqc_cmd_list(fastqc_path, fq, null_outdir), cmd_no_outdir_list)

    def test_form_fastqc_cmd_list_w_outdir_basic_example(self):
        """ test that form_fastq_cmd_list correctly generates a fastqc command list when passed valid arguments for
        the fastqc path, fastq input files, and output directory """

        # arguments to be formatted
        fastqc_path = '/usr/local/FastQC/./fastqc'
        fq = 'home/foo/bar/sample.fq'
        # option to redirect output not specified
        outdir = '/home/foo/bar/'

        cmd_w_outdir_list = ['/usr/local/FastQC/./fastqc', 'home/foo/bar/sample.fq', '--outdir', '/home/foo/bar/']

        self.assertEqual(exec_fastqc.form_fastqc_cmd_list(fastqc_path, fq, outdir), cmd_w_outdir_list)

    def test_form_fastqc_cmd_list_w_outdir_local_example(self):
        """ test that form_fastq_cmd_list correctly generates a fastqc command list when passed valid arguments for
        fastqc path, fastq input files, and output directory """

        # file and directory names to be used as output files
        sample_fq = 'sample_data.fastq'

        # join file paths to be relative to current directory
        sample_fq = os.path.join(self.input_dir, sample_fq)

        cmd_w_outdir_list = [self.fastqc_path, sample_fq, '--outdir', self.output_dir]

        self.assertEqual(exec_fastqc.form_fastqc_cmd_list(self.fastqc_path, sample_fq, self.output_dir), cmd_w_outdir_list)
    # endregion

    # region run_fastqc tests
    def test_run_fastqc_null_cmd_true(self):
        """ test that run_fastqc_cmd_list correctly raises a Value Error when invalid empty string(s) are used in place
        of required input """

        # arguments to be formatted
        null_fq_path = ''
        null_fq = ''
        # option for output directory omitted
        empty_outdir = None

        # catch a value exception when strings are empty
        with self.assertRaises(ValueError):
            exec_fastqc.run_fastqc(null_fq_path, null_fq, empty_outdir)

    def test_run_fastqc_short_input_good_stdout(self):
        """ test that fastqc subprocess call reports successful completion message to stdout when run_fastqc is
        passed valid arguments for the fastqc path, input fastq files, and output directory """

        # run fastqc and then check for expected stdout and stderr
        output, err = exec_fastqc.run_fastqc(self.fastqc_path, self.short_fq, self.output_dir)

        _, input_file = os.path.split(self.short_fq)

        fastqc_complete_str = "Analysis complete for {0}\n".format(input_file)
        # generate string expected of stdout
        self.assertEqual(fastqc_complete_str, output)

    def test_run_fastqc_output_zip_exists(self):
        """ test that a zipped directory contain fastqc reports exists when run_fastqc is passed valid arguments for
        the fastqc path, input fastq files, and output directory """

        # execute the function
        exec_fastqc.run_fastqc(self.fastqc_path, self.short_fq, self.output_dir)

        # obtain the input path
        _, input_file = os.path.split(self.short_fq)

        # trim off the file extension
        input_string_list = input_file.split(".")
        input_base = input_string_list[0]

        output_file = os.path.join(self.output_dir, input_base)

        fastqc_dir_zip = output_file + '_fastqc.zip'

        # test that the files exist are not empty
        self.assertTrue(os.stat(fastqc_dir_zip).st_size > 0)

    def test_run_fastqc_output_html_exists(self):
        """ test that an html report exists when run_fastqc is passed valid arguments for
        the fastqc path, input fastq files, and output directory """

        # execute the function
        exec_fastqc.run_fastqc(self.fastqc_path, self.short_fq, self.output_dir)

        # obtain the input path
        _, input_file = os.path.split(self.short_fq)

        # ignore the file extensions
        input_string_list = input_file.split(".")
        input_base = input_string_list[0]

        output_file = os.path.join(self.output_dir, input_base)

        fastqc_html = output_file + '_fastqc.html'

        # test that the files exist are not empty
        self.assertTrue(os.stat(fastqc_html).st_size > 0)
    #
    # endregion

    '''
    # TODO: Not concerned with the logging tests for now
    # TODO: Test passes when run individually but fails all tests are run in succession
    def test_run_fastqc_short_input_log_generated(self):
        logging.basicConfig(filename=self.log_file, level=logging.INFO)
        logging.debug('unittest: test_run_fastqc_short_input_log_generated')

        # run fastqc and then check for expected stdout and stderr
        output, err = exec_fastqc.run_fastqc(self.fastqc_path, self.short_fq, self.output_dir)

        print("normally i don't expect a null file: {0} ".format(self.log_file))

        # generate pipeline log
        self.assertTrue(os.stat(self.log_file).st_size > 0)
    '''

if __name__ == '__main__':
    unittest.main()
