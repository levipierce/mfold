from randomize_sequence import randomize_sequence as rs
import subprocess
import os
import shutil
import json
import sys 
import numpy as np
import argparse



def run(fname, seqs_to_gen, n_times_randomize_seq, mask=None):
    """
    mask should be a dictionary of keys with unigue names for the masked
    regions and relative positions of the regions i.e. 0 based
    """

    if mask:
        for key in mask:
            ps = process_sequence(fname, seqs_to_gen, n_times_randomize_seq, (key, mask[key]))
            save_results(fname, seqs_to_gen, ps, mask_name=key)
    else:
        ps = process_sequence(fname, seqs_to_gen, n_times_randomize_seq)
        save_results(fname, seqs_to_gen, ps)


def process_sequence(fname, n_seqs_to_generate, n_times_randomize_seq, mask=None):
    sample = {}
    for i in range(1, n_seqs_to_generate+1):
        fname_out = '%04d_%s' % (i, fname)
        seq_name = rs.process(fname, fname_out, n_times_randomize_seq, mask)
        base, ext = os.path.splitext(fname)
        if mask:
            fname_out = '%04d_%s_mask_%s%s' % (i, base, mask[0], ext)

        dir_name = "%04d" % i
        cmd = "python process_fasta.py --fasta_file %s --scratch_dir %s" % (fname_out, dir_name)
        print cmd
        #This gets rid of the vienna output
        with open(os.devnull, 'w') as FNULL:
            subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        #subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
        
        energy_file_name = "%s_energy.dat" % dir_name
        shutil.move("master_energy.dat", energy_file_name)
        #TODO Make this an option
        #shutil.rmtree(dir_name)
        with open(energy_file_name, "r") as fh:
            for l in fh:
                l_s = l.split()
                try:
                    sample[i].append((float(l_s[1]), float(l_s[2])))
                except KeyError:
                    sample[i] = [(float(l_s[1]), float(l_s[2]))]

    return sample

    
def save_results(fname, n_seqs_to_generate, sample, mask_name=None):
    """
    :param fname:
    :param n_seqs_to_generate:
    :param sample:
    :param mask_name:
    :return:
    """


    m_a = np.array([s[0] for s in sample[1]])
    m_b = np.array([s[1] for s in sample[1]])
    for i in range(2, n_seqs_to_generate+1):
        n_a = np.array([s[0] for s in sample[i]])
        n_b = np.array([s[1] for s in sample[i]])
        m_a = np.vstack((m_a, n_a))
        m_b = np.vstack((m_b, n_b))

    m_a_mean_of_segments = m_a.mean(0)
    m_b_mean_of_segments = m_b.mean(0)
    if n_seqs_to_generate != 1:
        m_a_mean_of_samples = m_a.mean(1)
        m_b_mean_of_samples = m_b.mean(1)

    fname_base = os.path.splitext(fname)[0]
    if mask_name:
        save_name = "%s_%s_mask_randomized_samples.json"%(fname_base, 
                                                            mask_name)
    else:
        save_name = "%s_randomized_samples.json"%(fname_base)

    with open(save_name, "w") as fh:
        json.dump(sample, fh)
        fh.write("\n")
        json.dump(m_a_mean_of_segments.tolist(), fh)
        fh.write("\n")
        json.dump(m_b_mean_of_segments.tolist(), fh)
        fh.write("\n")
        if n_seqs_to_generate != 1:
            json.dump(["%3f" % _ for _ in m_a_mean_of_samples.tolist()], fh)
            fh.write("\n")
            json.dump(["%3f" % _ for _ in m_b_mean_of_samples.tolist()], fh)
        else:
            json.dump(["0.00"], fh)
            fh.write("\n")
            json.dump(["0.00"], fh)

def read_mask(file_name):
    """
    :param file_name: name of the file containing the mask
    :return: dictionary of masked regions
    """
    '''
    Sample mask file...

    HSE1:		31795312-31795325
    HSE2:		31795405-31795418
    TATA:		31795486-31795491
    TSS: 		31795512-31795512
    TRIM28: 	31795571-31795578
    PB1: 		31795743-31795749
    PB2: 		31796358-31796364

    '''
    mask={}
    with open(file_name, "r") as fh:
        for l in fh:
            ls = l.split()
            name = ls[0].strip(":")
            s = int(ls[1].split("-")[0])
            e = int(ls[1].split("-")[1])
            try:
                mask[name].append((s, e))
            except KeyError:
                mask[name] = [(s, e)]
    return mask


def parse_it():
    parser = argparse.ArgumentParser(description='Fold randomized fasta file.')
    parser.add_argument(dest='fasta_file_name',
                        help='The name of the fasta file.')
    parser.add_argument('--num', dest='num_seqs_to_generate', type=int,
                        help='The number of randomized sequences to generate.', default=1)
    parser.add_argument('--num_rand', dest='num_rand',
                        help='The number of times to randomize a sequence.', default=1)
    parser.add_argument('--mask_file', dest='mask_file_name',
                        help='File containing regions masks.', default=None)
    parser.add_argument('--scratch_dir', dest='scratch_dir',
                        help='The name of the scratch dir', default=os.getcwd())
    return parser.parse_args()


if __name__=='__main__':
    args = parse_it()
    fasta_file_name = args.fasta_file_name
    seqs_to_gen = args.num_seqs_to_generate
    num_times_to_randomize = args.num_rand
    mask_file_name = args.mask_file_name


    run(fasta_file_name, seqs_to_gen, num_times_to_randomize, mask=mask_file_name)
    #run(fname_in, seqs_to_gen, num_times_to_randomize, 
    #    mask={"mask_test3":[(133450052,133450052+300), 
    #        (133450052+301,133450052+601)]})

