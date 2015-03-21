#!/usr/bin/python
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import argparse
import os
import copy

def read_sequence(fname):
    """
    """
    #Full path to file /Users/lpierce/mfold_development/CBFB.fasta 
    handle = open(fname, "rU")
    #store sequence
    sequence = []
    for record in SeqIO.parse(handle, "fasta") :
        #Build a list of short sequences:
        sequence = list(record.seq)
        header = record.description
        start_end_abs_pos = get_start_end_abs_position(header)

    handle.close()

    maps = get_position_maps(sequence, start_end_abs_pos)
    
    return (sequence, maps)

def get_start_end_abs_position(a_header):
    start, end = a_header.split(":")[1].split("-")
    return (int(start), int(end))

def get_position_maps(seq, start_end):
    rel_abs_pos_map = {}
    abs_rel_pos_map = {}
    start = start_end[0]

    for idx,n in enumerate(seq):
        rel_abs_pos_map[idx] = start + idx
        abs_rel_pos_map[start + idx] = idx

    return (rel_abs_pos_map, abs_rel_pos_map)


def shuffle_sequence(sequence, pos_maps, n=100,  mask=None):
    """
    chr6:31795202-31797300
    For supplied list of sites to mask need to make them relative to the
    range in the fasta
    """
    original_sequence_string = ''.join(sequence)
    temp_sequence = copy.copy(sequence)
    for i in range(n):
        random.shuffle(temp_sequence)
    #convert back to Seq object
    randomized_sequence_string = ''.join(temp_sequence)
    mask_name = "none"
    if mask:
        abs_rel_pos_map = pos_maps[1]
        mask_name = mask[0]
        #[(s0, e0), (s1,e1), (s2,e2)]
        mask_regions = mask[1]
        ran = list(randomized_sequence_string)
        ref = list(original_sequence_string)
        #can have multiple regions defined
        for region in mask_regions:
            s = abs_rel_pos_map[int(region[0])]
            e = abs_rel_pos_map[int(region[1]) + 1]
            for i in range(s,e):
                ran[i] = ref[i]
        randomized_sequence_string = ''.join(ran)
    #Double check first few nucleotides (The first 10 and last 10)
    print "Original sequence  : %s...%s " % (original_sequence_string[0:10], original_sequence_string[0:10])
    print "Radomized sequence : %s...%s " % (randomized_sequence_string[0:10], randomized_sequence_string[0:10])
    randomized_Seq_object = Seq(randomized_sequence_string, IUPAC.unambiguous_dna)

    randomized_sequence_record = SeqRecord(randomized_Seq_object,
                                           id="CBFB_randomized_mask_%s"%mask_name,
                                           name="CBFB",
                                           description="core-binding factor, beta subunit")
    return (randomized_sequence_record)


def write_new_sequence(fname, record, mask=None):
    """
    """
    #Write it out for mfold
    if mask:
        base, ext = os.path.splitext(fname)
        fname = "%s_mask_%s%s"%(base, mask[0], ext)
    output_handle = open(fname, "w")
    SeqIO.write(record, output_handle, "fasta")
    output_handle.close()
    print "Wrote randomized sequence to %s " % fname

def process(input_fname, output_fname, num, mask=None):
    """
    """
    seq, pos_maps = read_sequence(input_fname)
    randomized_record = shuffle_sequence(seq, pos_maps, num, mask)
    write_new_sequence(output_fname, randomized_record, mask)


def main(args):
    """
    """
    seq, header = read_sequence(args.input_fasta)
    randomized_record = shuffle_sequence(seq, args.num)
    write_new_sequence(args.output_fasta, randomized_record)
    

if __name__=="__main__":


    parser = argparse.ArgumentParser(description='Randomize fasta file')
    parser.add_argument('input_fasta', metavar='INPUT_FASTA_FILE', 
                   help='filename of the input fasta')
    parser.add_argument('output_fasta', metavar='OUTPUT_FASTA_FILE',
                   help='filename of the output fasta')
    parser.add_argument('--num', metavar='N', type=int, action='store',
                       default=100, help='number of times to randomize')

    args = parser.parse_args()
    print args
    main(args)
