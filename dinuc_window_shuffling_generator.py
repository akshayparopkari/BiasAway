# Luis del Peso
# Modified by A. Mathelier
# Vancouver, Jan 2012

import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from altschulEriksonDinuclShuffle import dinuclShuffle
from utils import GC, split_seq


def shuffle_window(ss, wl, step):
    bs = ss[:]
    for i in range(0, len(bs) - 1, step):
        # print i,"\t",ss[i:(i+wl)]
        bs = bs[0: i] + dinuclShuffle(bs[i: (i + wl)]) + bs[i + wl:]
    return(bs)   # returns shuffled sequence


def generate_sequences(seqs, winlen, step, nfold):
    bg_gc_list = []
    bg_lengths = []
    for record in seqs:
        seq = record.seq.__str__()
        for n in range(0, nfold):
            new_sequence = ""
            for sequence in split_seq(seq):
                if re.match("N", sequence):
                    new_sequence += sequence
                elif sequence:
                    new_sequence += shuffle_window(sequence, winlen, step)
            new_seq = SeqRecord(Seq(new_sequence, generic_dna),
                                id="background_seq_for_{}".format(record.name),
                                description="")
            print(new_seq.format("fasta"), end="")
            bg_gc_list.append(GC(new_sequence))
            bg_lengths.append(len(new_sequence))
    return bg_gc_list, bg_lengths
