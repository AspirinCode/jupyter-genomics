""" Module to call metaquast in pipeline"""

import argparse
import logging
import subprocess

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"

R1_MARKER = '_R1'
R2_MARKER = '_R2'


def form_metaquast_cmd_list(metaquast_fp, outdir, input_fasta):
    """format argument received to generate list to be used for metquast subprocess call

    Args:
        metaquast_fp(str): the string representing the path to metaquast executable
        outdir(str): the string representing the path to the output directory
        input_fasta(list): list of fasta files for the metaquast analysis

    Returns:
        call_args_list(list): is the sequence of commands to be passed to the metaquast call
    """

    if metaquast_fp is '' or metaquast_fp is None:
        raise ValueError('metaquast_path invalid. metaquast_path name is empty')
    if outdir is None:
        raise ValueError('outdir location invalid. outdir is None')
    if not input_fasta:
        raise ValueError('input contigs invalid. no fasta files specified')

    # required arguments
    call_args_list = ['python2', metaquast_fp]

    # add the fasta files
    call_args_list.extend(input_fasta)

    # add the output direcotry
    call_args_list.extend(['-o', outdir])

    return call_args_list


def run_metaquast(metaquast_fp, outdir, input_fasta):
    """ Call the metaquast software as a subprocess using input arguments

    Args:
        metaquast_fp(str): the string representing path to megahit executable
        outdir(str): the string representing the path to the output directory
        input_fq(list): list of input fasta files
    Returns:
        output(str): the string of characters printed by the metaquast call to stdout
        err(str): the string of characters rpinted by the metauqast call to stderr
    """

    logging.debug('beginning metaquast function call')

    # call the form_metaquast_cmd_list to generate the appropriate command
    call_args = form_metaquast_cmd_list(metaquast_fp, outdir, input_fasta)

    logging.info("calling popen with arguments {0}".format(" ".join(call_args)))

    # print(" ".join(call_args))
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def main():
    """ parse arguments and run_megahit using given string and provided fq files """

    # add basic arguments
    parser = argparse.ArgumentParser()
    # path to executable
    parser.add_argument("metaquast_path", type=str, help="path/to/metaquast_executable/")
    # select path to output directory, otherwise default location is input dir
    parser.add_argument("outdir", type=str, help="/path/to/output_directory")
    # select the contigs we wish to make the function call on
    parser.add_argument("input_fasta", nargs='+', type=str, help="input fasta contig files")

    # option to have log files
    parser.add_argument("-l", "--log", type=str, help="/path/to/log_file")

    args = parser.parse_args()

    # if a log file is specified, use it
    if args.log is not None:
        print("there exists an output log")
        log_file = args.log
        # begin program, set up logger
        logging.basicConfig(filename=log_file, level=logging.INFO)
        logging.debug('begin main')

    if args.log is not None:
        # sanity check
        logging.info('The metaquast executable file path is: {0}'.format(args.metaquast_path))
        logging.info('The pe2 fastq files processed by metaquast are: {0}'.format(args.input_fasta))
        logging.info('The megahit output directory is: {0}'.format(args.outdir))

    # call the megahit subprocess
    print(run_metaquast(args.metaquast_path, args.outdir, args.input_fasta))

if __name__ == "__main__":
    main()
