#!/usr/bin/python
#*-* coding: utf-8 *-*

###############################################################################
# BiasAway with the possibility of using very different ways of
# generating backgrounds lying into two categories:
# - Creation of new random sequences:
#   - di-nucleotide shuffling using the foreground sequences
#   - di-nucleotide shuffling within a sliding window using foreground
#   sequences
#- Extraction of sequences from a set of possible background sequences:
#   - respecting the %GC distribution of the foreground (using %GC bins)
#   - respecting the %GC distribution as in the previous item and also
#     respecting the %GC composition within a sliding window for %GC bin
###############################################################################

import argparse
import dinuc_shuffling_generator as dinuc_shuff
import dinuc_window_shuffling_generator as dinuc_win_shuff
import GC_compo_matching as GC_compo
import GC_window_compo_matching as GC_window_compo
from utils import get_seqs


def shuffling_generator(argu):
    seqs, __, __ = get_seqs(argu.fg_file)
    __, __ = dinuc_shuff.generate_sequences(seqs, argu.nfold)


def shuffling_window_generator(argu):
    seqs, __, __ = get_seqs(argu.fg_file)
    __, __ = dinuc_win_shuff.generate_sequences(seqs, argu.winlen,
                                                argu.step,
                                                argu.nfold)


def gc_compo_generator(argu):
    if argu.len_opt:
        gc_compo_len_generator(argu)
    else:
        gc_compo_generator_no_len(argu)


def gc_compo_generator_no_len(argu):
    __, fg_gc_bins, __ = GC_compo.fg_GC_bins(argu.fg_file)
    __, bg_gc_bins, __ = GC_compo.bg_GC_bins(argu.bg_file)
    __, __ = GC_compo.generate_sequences(fg_gc_bins, bg_gc_bins,
                                         argu.nfold)


def gc_compo_len_generator(argu):
    __, fg_gc_bins, __ = GC_compo.fg_len_GC_bins(argu.fg_file)
    __, bg_gc_bins, __ = GC_compo.bg_len_GC_bins(argu.bg_file)
    __, __ = GC_compo.generate_len_sequences(fg_gc_bins, bg_gc_bins,
                                             argu.nfold)


def gc_compo_window_generator(argu):
    if argu.len_opt:
        gc_compo_len_window_generator(argu)
    else:
        gc_compo_window_generator_no_len(argu)


def gc_compo_len_window_generator(argu):
    __, fg_gc_bins, __ = GC_window_compo.fg_len_GC_bins(argu.fg_file,
                                                        argu.winlen, argu.step)
    __, bg_gc_bins, __ = GC_window_compo.bg_len_GC_bins(argu.bg_file)
    __, __ = GC_window_compo.generate_len_sequences(fg_gc_bins, bg_gc_bins,
                                                    argu.deviation,
                                                    argu.winlen, argu.step,
                                                    argu.nfold)


def gc_compo_window_generator_no_len(argu):
    __, fg_gc_bins, __ = GC_window_compo.fg_GC_bins(argu.fg_file, argu.winlen,
                                                    argu.step)
    __, bg_gc_bins, __ = GC_window_compo.bg_GC_bins(argu.bg_file)
    __, __ = GC_window_compo.generate_sequences(fg_gc_bins, bg_gc_bins,
                                                argu.deviation, argu.winlen,
                                                argu.step, argu.nfold)


def shuffling_arg_parsing(subparsers):
    parser_d = subparsers.add_parser('d',
                                     help="di-nucleotide shuffling generator")
    parser_d.add_argument('-f', '--foreground', required=True, type=str,
                          dest="fg_file", action="store",
                          help="Foreground file in fasta format")
    help_str = "How many background sequences per each foreground sequence "
    help_str += "will be generated (default: 1)"
    parser_d.add_argument('-n', '--nfold', required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help=help_str)
    parser_d.set_defaults(func=shuffling_generator)


def window_shuffling_arg_parsing(subparsers):
    help_str = "di-nucleotide shuffling within a sliding window generator"
    parser_w = subparsers.add_parser('w', help=help_str)
    parser_w.add_argument("-w", "--winlen", required=False, type=int,
                          dest="winlen", action="store", default=200,
                          help="Window length (default: 200)")
    parser_w.add_argument("-s", "--step", required=False, type=int,
                          dest="step", action="store", default=1,
                          help="Sliding step (default: 1)")
    parser_w.add_argument('-f', '--foreground', required=True, type=str,
                          dest="fg_file", action="store",
                          help="Foreground file in fasta format")
    help_str = "How many background sequences per each foreground sequence "
    help_str += "will be generated (default: 1)"
    parser_w.add_argument('-n', '--nfold', required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help=help_str)
    parser_w.set_defaults(func=shuffling_window_generator)


def gc_compo_arg_parsing(subparsers):
    help_str = "%%GC distribution-based background chooser"
    parser_g = subparsers.add_parser('g', help=help_str)
    parser_g.add_argument("-b", "--background", required=True, type=str,
                          dest="bg_file", action="store",
                          help="Background file in fasta format")
    parser_g.add_argument('-f', '--foreground', required=True, type=str,
                          dest="fg_file", action="store",
                          help="Foreground file in fasta format")
    help_str = "How many background sequences per each foreground sequence "
    help_str += "will be choosen (default: 1)"
    parser_g.add_argument('-n', '--nfold', required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help=help_str)
    help_str = "Try to match the length as closely as possible "
    help_str += "(not set by default)"
    parser_g.add_argument('-l', '--length', required=False, dest="len_opt",
                          action="store_const", const=1, default=0,
                          help=help_str)
    parser_g.set_defaults(func=gc_compo_generator)


def gc_compo_window_arg_parsing(subparsers):
    help_str = "%%GC distribution and %%GC composition within a sliding "
    help_str += "window background chooser"
    parser_c = subparsers.add_parser('c', help=help_str)
    parser_c.add_argument("-b", "--background", required=True, type=str,
                          dest="bg_file", action="store",
                          help="Background file in fasta format")
    parser_c.add_argument("-w", "--winlen", required=False, type=int,
                          dest="winlen", action="store", default=200,
                          help="Window length (default: 200)")
    parser_c.add_argument("-s", "--step", required=False, type=int,
                          dest="step", action="store", default=1,
                          help="Sliding step (default: 1)")
    help_str = "Deviation from the mean (default: 2.6 for a "
    help_str += "threshold of mean + 2.6 * stdev)"
    parser_c.add_argument("-d", "--deviation", required=False, type=float,
                          dest="deviation", action="store", default=2.6,
                          help=help_str)
    parser_c.add_argument('-f', '--foreground', required=True, type=str,
                          dest="fg_file", action="store",
                          help="Foreground file in fasta format")
    help_str = "How many background sequences per each foreground sequence "
    help_str += "will be choosen (default: 1)"
    parser_c.add_argument('-n', '--nfold', required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help=help_str)
    help_str = "Try to match the length as closely as possible "
    help_str += "(not set by default)"
    parser_c.add_argument('-l', '--length', required=False, dest="len_opt",
                          action="store_const", const=1, default=0,
                          help=help_str)
    parser_c.set_defaults(func=gc_compo_window_generator)


def arg_parsing():
    descr = '''Background generator with the possibility of using very
    different ways of generating backgrounds lying into two categories:
        - Creation of new random sequences (generators):
            - di-nucleotide shuffling using the foreground sequences
            - di-nucleotide shuffling within a sliding window using foreground
                sequences
        - Extraction of sequences from a set of possible background sequences
          (choosers):
            - respecting the %GC distribution of the foreground (using %GC
              bins)
            - respecting the %GC distribution as in the previous item and also
            respecting the %GC composition within a sliding window for %GC bin
    '''
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(help="Choice of the generator/chooser",
                                       title="Subcommands",
                                       description="Valid subcommands")
    shuffling_arg_parsing(subparsers)
    window_shuffling_arg_parsing(subparsers)
    gc_compo_arg_parsing(subparsers)
    gc_compo_window_arg_parsing(subparsers)
    argu = parser.parse_args()
    return argu


###############################################################################
#                                   MAIN
###############################################################################
if __name__ == "__main__":
    arguments = arg_parsing()
    arguments.func(arguments)
