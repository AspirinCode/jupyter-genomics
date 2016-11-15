""" Module to call bowtie in pipeline """

import os
import subprocess
import argparse
import logging

__author__ = "YiDing Fang"
__maintainer__ = "YiDing Fang"
__email__ = "yif017@eng.ucsd.edu"
__status__ = "prototype"


def form_bowtie_build_cmd_list(bowtie_build_fp, input_contigs_fasta, output_index_fp):
    """ format arguments received to generate list used for bowtie_build subprocess call

    Args:
        bowtie_build_fp(str): the string representing the path to the bowtie program
        input_contigs_fasta(list): list of files which represent the fasta file inputs used for index construction
        output_index_fp(str): base name of file path to be used as output directory for index files

    Returns:
        call_args_list(list): the sequence of commands to be passed to the bowtie2 build call
    """
    if bowtie_build_fp is '':
        raise ValueError('bowtie2_build_path name is empty')
    if output_index_fp is '':
        raise ValueError('output file name invalid. index_output_fp is None')
    if input_contigs_fasta is '' or input_contigs_fasta is None:
        raise ValueError('no fasta file found')

    # required arguments
    calls_args_list = [bowtie_build_fp, input_contigs_fasta, output_index_fp]

    return calls_args_list


def form_bowtie_cmd_list(bowtie_fp, index_fp, pe1_fastq, pe2_fastq, u_fastq, output_sam_fp):
    """ format arguments recieved to generate list used for bowtie subprocess call

    Args:
        bowtie_fp(str): the string representing the path to the bowtie program
        index_fp(str): the string representing path to the bowtie index files
        pe1_fastq(list): list of files which represent the mate 1 forward reads
        pe2_fastq(list): list of files which represent the mate 2 reverse reads
        u_fastq(list): list of files containing unpaired reads to be aligned
        output_sam_fp(str):the string representing the path to the output sam file

    Returns:
        call_args(list): the sequence of commands to be passed to the bowtie2 call
    """

    # TODO this might not be necessary for now
    if bowtie_fp is '':
        raise ValueError('bowtie2_path is empty')
    if index_fp is '':
        raise ValueError('index_path is empty')
    if not pe1_fastq and not pe2_fastq and not u_fastq:
        raise ValueError('no fastq files specified')
    if output_sam_fp is '':
        raise ValueError('output_file_path is empty')

    # required arguments
    call_args_list = [bowtie_fp, '-x', index_fp]

    # add comma separated list of fastq files to process
    if pe1_fastq:
        call_args_list.append('-1')
        pe1_fq_str = ",".join(pe1_fastq)
        call_args_list.append(pe1_fq_str)
    if pe2_fastq:
        call_args_list.append('-2')
        pe2_fq_str = ",".join(pe2_fastq)
        call_args_list.append(pe2_fq_str)
    if u_fastq:
        call_args_list.append('-U')
        u_fq_str = ",".join(u_fastq)
        call_args_list.append(u_fq_str)

    call_args_list.extend(['-S', output_sam_fp])

    return call_args_list


def run_bowtie_build(bowtie_build_fp, input_contigs_fasta, output_index_fp):
    """Call the bowtie-build subprocess using input arguments to form index files

    Args:
        bowtie_build_fp(str): the string representing the path to the bowtie program
        input_contigs_fasta(list): list of files which represent the fasta file inputs and used for index construction
        output_index_fp(str): base name of file path to be used as output directory for index files

    Returns:
        output(str): the string of characters printed by the bowtie call to stdout
        err(str): the string of characters printed by the bowtie call to stderr
    """

    logging.debug('beginning bowtie-build function call')

    output_index_dir, base_name = os.path.split(output_index_fp)

    # call the form_bowtie_build_cmd_list method to generate the appropriate command
    call_args = form_bowtie_build_cmd_list(bowtie_build_fp, input_contigs_fasta, output_index_fp)

    if os.stat(output_index_dir):
        file_list = [f for f in os.listdir(output_index_dir) if f.endswith('.bt2')]
        if file_list:
            logging.warning("specified output directory for index contains bt2 files")
            raise OSError

    # write bowtie build command
    # print(call_args)
    logging.info("calling popen with arguments{0}".format(" ".join(call_args)))

    #print(" ".join(call_args)
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err


def run_bowtie(bowtie_fp, index_fp, pe1_fastq, pe2_fastq, u_fastq, output_sam_fp):
    """call the bowtie subprocess using input arguments to align reads with scaffold and generate output sam

    Args:
        bowtie_fp(str): the string representing the path to the bowtie program
        index_fp(str): the string representing path to the bowtie index files
        pe1_fastq(list): list of files which represent the mate 1 forward reads
        pe2_fastq(list): list of files which represent the mate 2 reverse reads
        u_fastq(list): list of files containing unpaired reads to be aligned
        output_sam_fp(str):the string representing the path to the output sam file

    Returns:
        output(str): the string of charcters printed by the bowtie call to stdout
        err(str): the string of characters rpinted by the bowtie call to stderr
    """

    logging.debug('beginning bowtie-align function call')

    # call the form_bowtie_cmd_list method to generate the appropriate command
    call_args = form_bowtie_cmd_list(bowtie_fp, index_fp, pe1_fastq, pe2_fastq, u_fastq, output_sam_fp)

    # write bowtie command
    # print(call_args)
    logging.info("calling popen with arguments{0}".format(" ".join(call_args)))

    #print(" ".join(call_args))
    process = subprocess.Popen(call_args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = process.communicate()

    return output, err

def build_run_bowtie(bowtie_fp, input_contigs_fasta, output_index_fp, pe1_fastq, pe2_fastq, u_fastq, output_sam_fp):
    """call run_bowtie_build and run_bowtie to generate output index from reference fasta and align fastq files to
    the generated contigs

    Args:
        bowtie_fp(str): the string representing the path to the bowtie program
        input_contigs_fasta(list): list of files which represent the fasta file inputs and used for index construction
        output_index_fp(str): base name of file path to be sued as output directory for index files
        pe1_fastq(list): list of files which represent the mate 1 forward reads
        pe2_fastq(list): list of files which represent the mate 2 reverse reads
        u_fastq(list): list of files containing unpaired reads to be aligned
        output_sam_fp(str):the( string representing the path to the output sam file

    Returns:
        buildout(str): the string of characters printed by the run_bowtie_build call to stdout
        builderr(str): the string of characters printed by the run_bowtie_build call to stderr

    """

    logging.debug('preparing bowtie build and align function calls')

    # bowtie build executable should be in the same directory as bowtie2 commands
    bowtie_build_fp = bowtie_fp + '-build'

    # build index files
    buildout, builderr = run_bowtie_build(bowtie_build_fp, input_contigs_fasta, output_index_fp)

    # once index files are built, run aligner using said index files
    stdout, stderr = run_bowtie(bowtie_fp, output_index_fp, pe1_fastq, pe2_fastq, u_fastq, output_sam_fp)

    return buildout, builderr, stdout, stderr


def main():
    """ parses arguments can calls run_megahit using given string and provided fq files """

    # add basic arguments
    parser = argparse.ArgumentParser()
    # path to executable
    parser.add_argument("bowtie_path", type=str, help="path/to/bowtie_executable/")
    # select path to output directory, otherwise default location is input dir
    parser.add_argument("index_output_fq", type=str, help="/path/to/output_index_directory/" )
    parser.add_argument("output_sam", type=str, help="path/to/output/sam_file")

    # select either paired end or single end
    parser.add_argument("-1", "--pe1", nargs='*', help="comma seperated list of forward pair fastq files")
    parser.add_argument("-2", "--pe2", nargs='*', help="comma seperated list of reverse pair fastq files")
    parser.add_argument("-u", "--unpaired", nargs='*', help="comma sperated list of unpaired fastq files")

    parser.add_argument("-c", "--contigs", nargs='+', required=True, help="comma seperated list of contig fasta files")

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
    logging.info('The bowtie executable file path is: {0}'.format(args.bowtie_path))

    if args.pe1:
        logging.info('The pe1 fastq files processed by bowtie are: {0}'.format(args.pe1))
    if args.pe2:
        logging.info('The pe2 fastq files processed by bowtie are: {0}'.format(args.pe2))
    if args.pe12:
        logging.info('The unpaired fastq files processed by bowtie are: {0}'.format(args.unpaired))
    if args.se:
        logging.info('The fasta files processed by bowtie-build are: {0}'.format(args.contigs))

    logging.info('The bowtie index directory is: {0}'.format(args.index_output_fq))
    logging.info('The bowtie output SAM file is: {0}'.format(args.output_sam))

    # call the megahit subprocess
    print(build_run_bowtie(args.bowtie_path, args.contigs, args.index_output_fq, args.pe1, args.pe2, args.unpaired,
                           args.output_sam))


if __name__ == "__main__":
    main()
