""" Module to call trimmomatic in pipeline"""

import os
import subprocess
import argparse
import logging
import re

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"

_R1_MARKER = '_R1'
_R2_MARKER = '_R2'
_ADAPTER_OPT_KNIGHT_LAB = 'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:70'


def form_trim_cmd_list(trimmomatic_fp, trim_mode, input_fastq, adapter_opt, outdir, trim_log):
    # TODO modify so that it can take either an adapter or a series of other trimmomatic options
    """ format arguments to be taken in by the fastqc function call. Handle logic to place output files in appropriate
        directories

    Args:
        trimmomatic_fp(str): the string representing the path to trimmomatic program
        trim_mode(str): the string informing us whether user is using paired or unpaired fastq files
        input_fastq(list): list of either 1 or 2 fq elements to be trimmed
        adapter_opt(
        outdir(str): the string representing the path to the output directory
        trim_log(str): option to include a trimlog

    Return value:
        call_args_list(list): is the sequence of commands to be passed to the trimmomatic call
    """

    if trimmomatic_fp is '':
        raise ValueError('trimmomatic_path name is empty')
    if trim_mode is not 'PE' and trim_mode is not 'SE':
        raise ValueError('unspecified if files are paired or unpaired')
    if adapter_opt is None or not adapter_opt:
        raise ValueError('adapter options not specified')
    if not input_fastq:
        raise ValueError('no fastq files specified')

    # required arguments
    call_args_list = ['java', '-jar', trimmomatic_fp, trim_mode]

    # add the fastq files to be trimmed
    if trim_mode is 'PE':
        if len(input_fastq) is 2:
            call_args_list.extend(input_fastq)
        else:
            raise ValueError('wrong number of fastq files for paired end reading')
    elif trim_mode is 'SE':
        if len(input_fastq) is 1:
            call_args_list.extend(input_fastq)
        else:
            raise ValueError('wrong number of fastq files for single end reading')

    # default location for trimmed files is in current directory if not split
    indir, input_file = os.path.split(input_fastq[0])

    # remove paired ending extensions from fastq file to form base output name
    if _R1_MARKER in input_file:
        base_fq = re.sub(_R1_MARKER, '', input_file, 1)
    elif _R2_MARKER in input_file:
        base_fq = re.sub(_R2_MARKER, '', input_file, 1)
    else:
        base_fq = input_file

    # extend output file names to complete paths
    if outdir is not None:
        base_out = os.path.join(outdir, base_fq)
    else:
        base_out = os.path.join(indir, base_fq)

    # output placed in specified outdir
    if trim_mode is 'PE':
        call_args_list.extend(["-baseout", base_out])
    elif trim_mode is 'SE':
        call_args_list.extend([base_out])

    # add trimlog option if selected
    if trim_log is not None:
        call_args_list.extend(['-trimlog', trim_log])

    # add fastqc usage settings
    call_args_list.extend(adapter_opt)

    return call_args_list


def run_trimmomatic(trimmomatic_fp, trim_mode, input_fastq, adapter_opt, outdir, trim_log):
    """ Call the fastqc software as a subprocess using input arguments

    Args:
        trimmomatic_fp(str): the string representing path to fastqc program
        trim_mode(str): either PE or SE depending on user inputs
        input_fastq(list): list containing the potential fastq files to be used
        adapter_opt(list): formatted sequence for settings
        outdir(str): the string representing the path to the output directory
        trim_log(str): name of the file to write trimmomatic log to

    Return value:
        output(str): the string of characters printed by the trimmomatic call to std out
        err(str): the string of characters printed by the trimmomatic call to std err
    """

    logging.debug('beginning run_trimmomatic function call')
    # call the form_trimmomatic_cmd_list method to generate the appropriate command
    call_args = form_trim_cmd_list(trimmomatic_fp, trim_mode, input_fastq, adapter_opt, outdir, trim_log)

    # write trimmomatic command
    logging.info("calling popen with arguments {0}".format(" ".join(call_args)))

    # print(" ".join(call_args))
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def run_pipeline_trimmomatic(trimmomatic_fp, input_fastq, outdir, trim_log):
    """ Call the fastqc software as a subprocess using input arguments and presets given for this pipeline

    Args:
        trimmomatic_fp(str): the string representing path to fastqc program
        input_fastq(list): list containing the potential fastq files to be used
        outdir(str): the string representing the path to the output directory
        trim_log(str): name of the file to write trimmomatic log to

    Return value:
        output(str): the string of characters printed by the trimmomatic call to std out
        err(str): the string of characters printed by the trimmomatic call to std err
    """

    # check if a log needs to be written
    logging.debug('beginning run_trimmomatic function call')
    # call the form fastqc line method to generate the appropriate command
    call_args = form_trim_cmd_list(trimmomatic_fp, 'PE', input_fastq, _ADAPTER_OPT_KNIGHT_LAB, outdir, trim_log)

    # write trimmomatic command
    logging.info("calling popen with arguments {0}".format(" ".join(call_args)))

    # print(" ".join(call_args))
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def main():
    """ parse arguments and calls run_trimmomatic using given settings and fq files """

    # add basic arguments
    trim_mode = ''
    input_fq = []

    parser = argparse.ArgumentParser()
    # path to executable
    parser.add_argument("trimmomatic_path", type=str, help="path/to/trimmomatic.jar/")
    # select either paired end or single end
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-PE", "--paired_end", type=str, nargs=2, help="/path/to/paired_fastq_files")
    group.add_argument("-SE", "--single_end", type=str, nargs=1, help="/path/to/unpaired_fastq_files")
    # must input valid adapter options
    parser.add_argument("adapter_opt", type=str, nargs='+', help="trim settings as defined by trimmomatic")
    # select path to output directory, otherwise default location is input dir
    parser.add_argument("-o", "--outdir", type=str, help="/path/to/output_directory")
    # option to have log files
    parser.add_argument("-l", "--log", type=str, help="/path/to/log_file")
    parser.add_argument("-L", "--trim_log", type=str, help="/path/to/trim_log")

    args = parser.parse_args()

    if args.paired_end:
        trim_mode = 'PE'
        input_fq = args.paired_end
    if args.single_end:
        trim_mode = 'SE'
        input_fq = args.single_end

    # if a log file is specified, use it
    if args.log is not None:
        print("there exists an output log")
        log_file = args.log
        # begin program, set up logger
        logging.basicConfig(filename=log_file, level=logging.INFO)

    logging.debug('begin main')

    logging.info('The trimmomatic jar file path is: {0}'.format(args.trimmomatic_path))

    if args.paired_end:
        print('The paired fastq files processed by trimmomatic is: {0}, {1}'.format(args.paired_end[0],
                                                                                    args.paired_end[1]))
        logging.info('The paired fastq files processed by trimmomatic is: {0}, {1}'.format(args.paired_end[0],
                                                                                           args.paired_end[1]))
    elif args.single_end:
        print('The unpaired fastq file processed by trimmomatic is: {0}'.format(args.paired_end[0]))
        logging.info('The paired fastq files processed by trimmomatic is: {0}'.format(args.paired_end[0]))

    logging.info('The trimmomatic settings are {0}'.format(args.adapter_opt))
    logging.info('The trimmomatic output directory is: {0}'.format(args.outdir))
    logging.info('The trim log file is: {0}'.format(args.trim_log))

    # call the trimmomatic subprocess
    print(run_trimmomatic(args.trimmomatic_path, trim_mode, input_fq, args.adapter_opt, args.outdir, args.trim_log))


if __name__ == "__main__":
    main()