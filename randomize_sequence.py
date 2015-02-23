
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import argparse

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
    handle.close()

    
    return sequence


def shuffle_sequence(sequence, n=100):
    """
    """
    original_sequence_string = ''.join(sequence)
    temp_sequence = sequence
    for i in range(100):
        random.shuffle(temp_sequence)
    #convert back to Seq object
    randomized_sequence_string = ''.join(temp_sequence)
    #Double check first few nucleotides (The first 10 and last 10)
    print "Original sequence  : %s...%s " % (original_sequence_string[0:10], original_sequence_string[0:10])
    print "Radomized sequence : %s...%s " % (randomized_sequence_string[0:10], randomized_sequence_string[0:10])
    randomized_Seq_object = Seq(randomized_sequence_string, IUPAC.unambiguous_dna)

    randomized_sequence_record = SeqRecord(randomized_Seq_object,
                                           id="CBFB_randomized",
                                           name="CBFB",
                                           description="core-binding factor, beta subunit")
    return randomized_sequence_record


def write_new_sequence(fname, record):
    """
    """
    #Write it out for mfold
    output_handle = open(fname, "w")
    SeqIO.write(record, output_handle, "fasta")
    output_handle.close()
    print "Wrote randomized sequence to %s " % fname

def main(args):
    """
    """
    seq = read_sequence(args.input_fasta)
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
