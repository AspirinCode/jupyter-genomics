""" Module to call fastqc in pipeline"""

import os
import subprocess
import argparse
import logging

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"


def form_fastqc_cmd_list(fastqc_fp, fastq_fp, outdir):
    """Generate argument list to be given as input to the fastqc function call.

    Args:
        fastqc_fp(str): the string representing path to fastqc program
        fastq_fp(str): the string representing path to the fastq file to be evaluated
        outdir(str): the string representing the path to the output directory

    Return value:
        call_args(list): the list of call_args representing the options for the fastqc subprocess call

    Raises:
        ValueError is raised when either the fastqc path or the fastqc input files are empty
    """
    # throw exceptions to prevent user from accidentally using interactive fastqc
    if fastqc_fp is '':
        raise ValueError('fastqc_fp name is empty')
    if fastq_fp is '':
        raise ValueError('fastq_fp file name is empty')

    # required arguments
    call_args_list = [fastqc_fp, fastq_fp]

    # direct output
    if outdir is not None:
        call_args_list.extend(["--outdir", outdir])

    return call_args_list


def run_fastqc(fastqc_fp, fastq_fp, outdir):
    """ Call the fastqc software as a subprocess using input arguments.

    Args:
        fastqc_fp(str): the string representing path to fastqc program
        fastq_fp(str): the string representing the file to be evaluated
        outdir(str): the string representing the path to the output directory

    Return value:
        output(str): the string of characters sent to stdout by fastqc
        err(str): the string of characters sent to stderr by fastq
    """
    logging.debug('beginning run_fastqc function call')

    # call the form_fastqc_cmd_list method to generate the appropriate command
    call_args = form_fastqc_cmd_list(fastqc_fp, fastq_fp, outdir)
    logging.info("calling popen with arguments '{0}'".format(" ".join(call_args)))

    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def main():
    """Parse command-line arguments and pass them to the fastqc software."""

    parser = argparse.ArgumentParser()
    parser.add_argument("fastqc", type=str, help="path/to/fastqc_executable/")
    parser.add_argument("fastq_fp", type=str, help="path/to/fastq")
    parser.add_argument("-o", "--outdir", type=str, help="/path/to/output_directory")
    parser.add_argument("-l", "--log", type=str, help="/path/to/log_file")
    args = parser.parse_args()

    if args.outdir is not None:
        output_dir = args.outdir
    else:
        # use the directory of the input file for output
        output_dir, _ = os.path.split(args.fastq_fp)

    # if a log file is specified, set up the logger
    if args.log is not None:
        log_file = args.log
        logging.basicConfig(filename=log_file, level=logging.INFO)

    logging.debug('begin main')
    logging.info('The fastqc file path is: %s' % args.fastqc)
    logging.info('The fastq file processed by fastqc is: %s' % args.fastq_fp)
    logging.info('The fastqc output directory is: %s' % args.outdir)

    # call the fastqc subprocess
    print(run_fastqc(args.fastqc, args.fastq_fp, output_dir))


if __name__ == "__main__":
    main()
