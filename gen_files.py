from randomize_sequence import process
import subprocess
import os
import shutil
import json
import sys 
import numpy as np


def run(fname, seqs_to_gen, n_times_randomize_seq, mask=None):
    """
    mask should be a dictionary of keys with unigue names for the masked
    regions and relative positions of the regions i.e. 0 based
    """

    if mask:
        for key in mask:
            ps = process_sequence(fname, seqs_to_gen, n_times_randomize_seq, (key,mask[key]))
            save_results(fname, seqs_to_gen, ps, mask_name=key)
    else:
        ps = process_sequence(fname, seqs_to_gen, n_times_randomize_seq)
        save_results(fname, seqs_to_gen, ps)


def process_sequence(fname, n_seqs_to_generate, n_times_randomize_seq, mask=None):
    sample = {}
    for i in range(1, n_seqs_to_generate+1):

        fname_out ='%04d_%s' % (i, fname)
        process(fname, fname_out, n_times_randomize_seq, mask)
        fname_in = fname_out

        dir_name = "%04d"%i
        cmd = "python process_fasta.py --fasta_file %s --scratch_dir blah/%s"% (fname_out, dir_name)
        print cmd
        with open(os.devnull, 'w') as FNULL:
            subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        
        energy_file_name = "%s_energy.dat"%dir_name
        shutil.move("master_energy.dat", energy_file_name)
        shutil.rmtree(dir_name)
        with open(energy_file_name, "r") as fh:
            for l in fh:
                l_s = l.split()
                try:
                    sample[i].append((float(l_s[1]),float(l_s[2])))
                except KeyError:
                    sample[i] = [(float(l_s[1]),float(l_s[2]))]

    return sample

    
def save_results(fname, n_seqs_to_generate, sample, mask_name=None):


    m_a = np.array([s[0] for s in sample[1]])
    m_b = np.array([s[1] for s in sample[1]])
    for i in range(2, n_seqs_to_generate+1):
        n_a = np.array([s[0] for s in sample[i]])
        n_b = np.array([s[1] for s in sample[i]])
        m_a = np.vstack((m_a,n_a))
        m_b = np.vstack((m_b,n_b))

    m_a_mean_of_segments = m_a.mean(0)
    m_b_mean_of_segments = m_b.mean(0)
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
        json.dump(m_a_mean_of_segments.tolist(),fh)
        fh.write("\n")
        json.dump(m_b_mean_of_segments.tolist(),fh)
        fh.write("\n")
        json.dump(["%3f"%_ for _ in m_a_mean_of_samples.tolist()],fh)
        fh.write("\n")
        json.dump(m_b_mean_of_samples.tolist(),fh)



if __name__=='__main__':
    fname_in = sys.argv[1]
    seqs_to_gen = int(sys.argv[2])
    num_times_to_randomize = int(sys.argv[3])

    run(fname_in, seqs_to_gen, num_times_to_randomize, 
        mask={"mask_test3":[(133450052,133450052+300), (133450052+301,133450052+601)]})
