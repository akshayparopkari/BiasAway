#!/usr/bin/python

"""
BiasAway module generating adapted background for motif overrepresentation.

 BiasAway with the possibility of using very different ways of
 generating backgrounds lying into two categories:
 - Creation of new random sequences:
   - mono-nucleotide shuffling using the foreground sequences
   - mono-nucleotide shuffling within a sliding window using foreground
     sequences
   - di-nucleotide shuffling using the foreground sequences
   - di-nucleotide shuffling within a sliding window using foreground
     sequences
- Extraction of sequences from a set of possible background sequences:
   - respecting the GC distribution of the foreground (using GC bins)
   - respecting the GC distribution as in the previous item and also
     respecting the GC composition within a sliding window for GC bin

"""

import argparse
import mononuc_shuffling_generator as mononuc_shuff
import mononuc_window_shuffling_generator as mononuc_win_shuff
import dinuc_shuffling_generator as dinuc_shuff
import dinuc_window_shuffling_generator as dinuc_win_shuff
import GC_compo_matching as GC_compo
import GC_window_compo_matching as GC_window_compo
from utils import get_seqs
import sys
import os
import errno


def mononuc_shuffling_generator(argu):
    seqs, _, _ = get_seqs(argu.fg_file)
    _, _ = mononuc_shuff.generate_sequences(seqs, argu.nfold)


def dinuc_shuffling_generator(argu):
    seqs, _, _ = get_seqs(argu.fg_file)
    _, _ = dinuc_shuff.generate_sequences(seqs, argu.nfold)


def mononuc_shuffling_window_generator(argu):
    seqs, _, _ = get_seqs(argu.fg_file)
    _, _ = mononuc_win_shuff.generate_sequences(seqs, argu.winlen, argu.step, argu.nfold)


def dinuc_shuffling_window_generator(argu):
    seqs, _, _ = get_seqs(argu.fg_file)
    _, _ = dinuc_win_shuff.generate_sequences(seqs, argu.winlen, argu.step, argu.nfold)


def test_empty_bg_dir(bg_dir):
    if os.path.isdir(bg_dir):
        if os.listdir(bg_dir):
            sys.exit("EXITING since both a non-empty background directory and a "
                     "background file are given")
    else:
        try:
            os.makedirs(bg_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(bg_dir):
                pass
            else:
                raise


def test_non_empty_bg_dir(bg_dir):
    if not (os.path.isdir(bg_dir) and os.listdir(bg_dir)):
        sys.exit("EXITING since the background directory does not exist or is empty")


def gc_compo_generator(argu):
    if argu.len_opt:
        gc_compo_len_generator(argu)
    else:
        gc_compo_generator_no_len(argu)


def gc_compo_generator_no_len(argu):
    _, fg_gc_bins, _ = GC_compo.fg_GC_bins(argu.fg_file)
    bg_gc_bins = None
    if argu.bg_file:
        test_empty_bg_dir(argu.bg_dir)
        _, bg_gc_bins, _ = GC_compo.bg_GC_bins(argu.bg_file, argu.bg_dir)
    else:
        test_non_empty_bg_dir(argu.bg_dir)
    _, _ = GC_compo.generate_sequences(fg_gc_bins, bg_gc_bins, argu.bg_dir, argu.nfold)


def gc_compo_len_generator(argu):
    _, fg_gc_bins, _ = GC_compo.fg_len_GC_bins(argu.fg_file)
    bg_gc_bins = None
    if argu.bg_file:
        test_empty_bg_dir(argu.bg_dir)
        _, bg_gc_bins, _ = GC_compo.bg_len_GC_bins(argu.bg_file, argu.bg_dir)
    else:
        test_non_empty_bg_dir(argu.bg_dir)
    _, _ = GC_compo.generate_len_sequences(fg_gc_bins, bg_gc_bins, argu.bg_dir,
                                           argu.nfold)


def gc_compo_window_generator(argu):
    if argu.len_opt:
        gc_compo_len_window_generator(argu)
    else:
        gc_compo_window_generator_no_len(argu)


def gc_compo_len_window_generator(argu):
    _, fg_gc_bins, _ = GC_window_compo.fg_len_GC_bins(argu.fg_file, argu.winlen,
                                                      argu.step)
    bg_gc_bins = None
    if argu.bg_file:
        test_empty_bg_dir(argu.bg_dir)
        _, bg_gc_bins, _ = GC_window_compo.bg_len_GC_bins(argu.bg_file, argu.bg_dir)
    else:
        test_non_empty_bg_dir(argu.bg_dir)
    _, _ = GC_window_compo.generate_len_sequences(fg_gc_bins, bg_gc_bins, argu.bg_dir,
                                                  argu.deviation, argu.winlen, argu.step,
                                                  argu.nfold)


def gc_compo_window_generator_no_len(argu):
    _, fg_gc_bins, _ = GC_window_compo.fg_GC_bins(argu.fg_file, argu.winlen, argu.step)
    bg_gc_bins = None
    if argu.bg_file:
        test_empty_bg_dir(argu.bg_dir)
        _, bg_gc_bins, _ = GC_window_compo.bg_GC_bins(argu.bg_file, argu.bg_dir)
    else:
        test_non_empty_bg_dir(argu.bg_dir)
    _, _ = GC_window_compo.generate_sequences(fg_gc_bins, bg_gc_bins, argu.bg_dir,
                                              argu.deviation, argu.winlen, argu.step,
                                              argu.nfold)


def mononuc_shuffling_arg_parsing(subparsers):
    parser_d = subparsers.add_parser("m", help="mono-nucleotide shuffling generator")
    parser_d.add_argument("-f", "--foreground", required=True, type=str, dest="fg_file",
                          action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_d.add_argument("-n", "--nfold", required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be generated (default: 1)")
    parser_d.set_defaults(func=mononuc_shuffling_generator)


def mononuc_window_shuffling_arg_parsing(subparsers):
    parser_w = subparsers.add_parser("f", help="mono-nucleotide shuffling within a "
                                     "sliding window generator")
    parser_w.add_argument("-w", "--winlen", required=False, type=int, dest="winlen",
                          action="store", default=100,
                          help="Window length (default: 100 bp)")
    parser_w.add_argument("-s", "--step", required=False, type=int, dest="step",
                          action="store", default=1,  help="Sliding step (default: 1 bp)")
    parser_w.add_argument("-f", "--foreground", required=True, type=str, dest="fg_file",
                          action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_w.add_argument("-n", "--nfold", required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be generated (default: 1 bp)")
    parser_w.set_defaults(func=mononuc_shuffling_window_generator)


def dinuc_shuffling_arg_parsing(subparsers):
    parser_d = subparsers.add_parser("d", help="di-nucleotide shuffling generator")
    parser_d.add_argument("-f", "--foreground", required=True, type=str, dest="fg_file",
                          action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_d.add_argument("-n", "--nfold", required=False, type=int, dest="nfold",
                          action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be generated (default: 1)")
    parser_d.set_defaults(func=dinuc_shuffling_generator)


def dinuc_window_shuffling_arg_parsing(subparsers):
    parser_w = subparsers.add_parser("w", help="di-nucleotide shuffling within a sliding "
                                     "window generator")
    parser_w.add_argument("-w", "--winlen", required=False, type=int, dest="winlen",
                          action="store", default=100,
                          help="Window length (default: 100 bp)")
    parser_w.add_argument("-s", "--step", required=False, type=int, dest="step",
                          action="store", default=1, help="Sliding step (default: 1bp)")
    parser_w.add_argument("-f", "--foreground", required=True, type=str, dest="fg_file",
                          action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_w.add_argument("-n", "--nfold", required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be generated (default: 1 bp)")
    parser_w.set_defaults(func=dinuc_shuffling_window_generator)


def gc_compo_arg_parsing(subparsers):
    parser_g = subparsers.add_parser("g", help="GC content-based background chooser")
    parser_g.add_argument("-r", "--bgdirectory", required=True, type=str,
                          dest="bg_dir", action="store",
                          help="Background directory [REQUIRED]")
    parser_g.add_argument("-b", "--background", required=False, type=str,
                          dest="bg_file", action="store",
                          help="Background file in FASTA format")
    parser_g.add_argument("-f", "--foreground", required=True, type=str,
                          dest="fg_file", action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_g.add_argument("-n", "--nfold", required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be choosen (default: 1 bp)")
    parser_g.add_argument("-l", "--length", required=False, dest="len_opt",
                          action="store_const", const=1, default=0,
                          help="Try to match the length as closely as possible")
    parser_g.set_defaults(func=gc_compo_generator)


def gc_compo_window_arg_parsing(subparsers):
    parser_c = subparsers.add_parser("c", help="GC distribution and GC content within"
                                     "a sliding window background chooser")
    parser_c.add_argument("-r", "--bgdirectory", required=True, type=str, dest="bg_dir",
                          action="store", help="Background directory [REQUIRED]")
    parser_c.add_argument("-b", "--background", required=False, type=str, dest="bg_file",
                          action="store", help="Background file in FASTA format")
    parser_c.add_argument("-w", "--winlen", required=False, type=int, dest="winlen",
                          action="store", default=100,
                          help="Window length (default: 100 bp)")
    parser_c.add_argument("-s", "--step", required=False, type=int, dest="step",
                          action="store", default=1, help="Sliding step (default: 1 bp)")
    parser_c.add_argument("-d", "--deviation", required=False, type=float,
                          dest="deviation", action="store", default=2.6,
                          help="Deviation from the mean (default: 2.6 for a threshold of "
                          "mean + 2.6 * stdev)")
    parser_c.add_argument("-f", "--foreground", required=True, type=str, dest="fg_file",
                          action="store",
                          help="Foreground file in FASTA format [REQUIRED]")
    parser_c.add_argument("-n", "--nfold", required=False, type=int,
                          dest="nfold", action="store", default=1,
                          help="How many background sequences per each foreground "
                          "sequence will be choosen (default: 1 bp)")
    parser_c.add_argument("-l", "--length", required=False, dest="len_opt",
                          action="store_const", const=1, default=0,
                          help="Try to match the length as closely as possible "
                          "(not set by default)")
    parser_c.set_defaults(func=gc_compo_window_generator)


def arg_parsing():
    descr = """Background generator with the possibility of using very different ways of
    generating backgrounds lying into two categories:
        - Creation of new random sequences (generators):
            - mono-nucleotide shuffling using the foreground sequences
            - mono-nucleotide shuffling within a sliding window using foreground sequences
            - di-nucleotide shuffling using the foreground sequences
            - di-nucleotide shuffling within a sliding window using foreground sequences
        - Extraction of sequences from a set of possible background sequences (choosers):
            - respecting the GC distribution of the foreground (using GC bins)
            - respecting the GC distribution as in the previous item and also respecting
              the GC composition within a sliding window for GC bin
    """
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(help="Choice of the generator/chooser",
                                       title="Subcommands",
                                       description="Valid subcommands")
    mononuc_shuffling_arg_parsing(subparsers)
    mononuc_window_shuffling_arg_parsing(subparsers)
    dinuc_shuffling_arg_parsing(subparsers)
    dinuc_window_shuffling_arg_parsing(subparsers)
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
