"""
Generation of sequences using mono-nucleotide of input sequences.

Module allowing the generation of sequences by using a mono-nucleotide.
shuffling of the given sequences.

"""

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from utils import GC
import random


def generate_sequences(seqs, nfold):
    """
    Generate sequences by shuffling input (mononucleotide).

    return tuple containing %GC compo and length distrib of output.

    """
    bg_gc_list = []
    bg_lengths = []
    for record in seqs:
        seq = record.seq.__str__()
        for _ in range(0, nfold):
            new_sequence = "".join(random.sample(seq, len(seq)))
            new_seq = SeqRecord(Seq(new_sequence, generic_dna),
                                id="background_seq_for_{}".format(record.name),
                                description="")
            print(new_seq.format("fasta"), end="")
            bg_gc_list.append(GC(new_sequence))
            bg_lengths.append(len(new_sequence))
    return bg_gc_list, bg_lengths
