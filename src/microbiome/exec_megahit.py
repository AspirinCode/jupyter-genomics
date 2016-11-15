""" Module to call megahit in pipeline """

import argparse
import logging
import subprocess

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"


def form_megahit_cmd_list(megahit_fp, outdir, max_memory, pe1_fastq, pe2_fastq, pe12_fastq, se_fastq):
    """ format arguments received to generate list to be used for megahit subprocess call

    Args:
        megahit_fp(str): the string representing path to megahit executable
        max_memory(str): maximum bytes of memory to be used for construction of sparse de Bruijn graph
        pe1_fastq(list): list of paired foward fastq files
        pe2_fastq(list): list of paired reverse fastq files
        pe12_fastq(list): list of interleaved fastq files
        se_fastq(list): list of single ended fastq files
        outdir(str): a string representing the path to the output directory

    Returns:
        call_args_list(list): is the sequence of commands to be passed to the megahit call
    """

    # TODO this might not be necessary for now
    if megahit_fp is '':
        raise ValueError('megahit_path name is empty')
    if outdir is None:
        raise ValueError('outdir location invalid. outdir is None')
    if not pe1_fastq and not pe2_fastq and not pe12_fastq and not se_fastq:
        raise ValueError('no fastq files specified')

    # required arguments
    call_args_list = [megahit_fp]

    # add fastqc usage settings
    call_args_list.extend(['-o', outdir])

    if max_memory is not None:
        call_args_list.extend(['-m', max_memory])

    # add comma separated list of fastq files to process
    if pe1_fastq:
        call_args_list.append('-1')
        pe1_fq_str = ",".join(pe1_fastq)
        call_args_list.append(pe1_fq_str)
    if pe2_fastq:
        call_args_list.append('-2')
        pe2_fq_str = ",".join(pe2_fastq)
        call_args_list.append(pe2_fq_str)
    if pe12_fastq:
        call_args_list.append('--12')
        pe12_fq_str = ",".join(pe12_fastq)
        call_args_list.append(pe12_fq_str)
    if se_fastq:
        call_args_list.append('-r')
        se_fq_str = ",".join(se_fastq)
        call_args_list.append(se_fq_str)

    return call_args_list


def run_megahit(megahit_fp, outdir, max_memory, pe1_fastq, pe2_fastq, pe12_fastq, se_fastq):
    """ Call the megahit software subproccess using input arguments

    Args:
        megahit_fp(str): the string representing path to megahit executable
        max_memory(str): maximum bytes of memory to be used for construction of sparse de Bruijn graph
        pe1_fastq(list): list of paired foward fastq files
        pe2_fastq(list): list of paired reverse fastq files
        pe12_fastq(list): list of interleaved fastq files
        se_fastq(list): list of single ended fastq files
        outdir(str): a string representing the path to the output directory

    Returns:
        output(str): the string of characters printed by the megahit call to std out
        err(str): the string of characters printed by the megahit call to std err
    """

    logging.debug('beginning megahit function call')

    # call the form_megahit_cmd_list line method to generate the appropriate command
    call_args = form_megahit_cmd_list(megahit_fp, outdir, max_memory, pe1_fastq, pe2_fastq, pe12_fastq,
                                      se_fastq)

    # write megahit command
    # print(call_args)
    logging.info("calling popen with arguments {0}".format(" ".join(call_args)))

    # print(" ".join(call_args))
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def main():
    """
    parse arguments and call run_megahit using given string and provided fq files
    """
    # add basic arguments
    parser = argparse.ArgumentParser()
    # path to executable
    parser.add_argument("megahit_path", type=str, help="path/to/megahit_executable/")
    # select path to output directory, otherwise default location is input dir
    parser.add_argument("outdir", type=str, help="/path/to/output_directory")

    # memory option
    parser.add_argument("-m", "--memory", type=str, help="maximum bytes of memory to be used in graph construction")

    # select either paired end or single end
    parser.add_argument("-1", "--pe1", nargs='*', help="forward pair fastq")
    parser.add_argument("-2", "--pe2", nargs='*', help="reverse pair fastq")
    parser.add_argument("-12", "--pe12", nargs='*', help="interleaved fastq")
    parser.add_argument("-r", "--se", nargs='*', help="single-end fastq")

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

    # sanity check
    logging.info('The megahit executable file path is: {0}'.format(args.megahit_path))

    if args.pe1:
        logging.info('The pe1 fastq files processed by megahit are: {0}'.format(args.pe1))
    if args.pe2:
        logging.info('The pe2 fastq files processed by megahit are: {0}'.format(args.pe2))
    if args.pe12:
        logging.info('The pe12 fastq files processed by megahit are: {0}'.format(args.pe12))
    if args.se:
        logging.info('The pe2 fastq files processed by megahit are: {0}'.format(args.se))

    logging.info('The megahit output directory is: {0}'.format(args.outdir))

    # call the megahit subprocess
    print(run_megahit(args.megahit_path, args.outdir, args.memory, args.pe1, args.pe2, args.pe12, args.se))


if __name__ == "__main__":
    main()
