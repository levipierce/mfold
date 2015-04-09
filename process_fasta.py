#!/usr/bin/python
import tarfile
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
import shutil
import argparse
import time
import logging
import glob

os.environ["VIENNA"] = "/home/lpierce/Software/vienna/share/ViennaRNA"
os.environ["VIENNA"] = "/home/wanglab/software/viennaRNA-2.1.5/share/ViennaRNA"
#These should be settings initalized when creating the class
FOLD_SIZE = 300
SLIDING_WINDOW_SIZE = 150
NUM_PROCESSORS = 4

def parse_it():
    parser = argparse.ArgumentParser(description='Process a fasta file.')
    parser.add_argument('--fasta_file', dest='fasta_file_name',help='The name of the fasta file.')
    parser.add_argument('--scratch_dir', dest='scratch_dir', help='The name of the scratch dir')
    args = parser.parse_args()
    return args

def run_vienna(sequence):                                                       
    """                                                                         
    """                                                                         
    #RNAfold  --gquad   --noconv   <    test_fasta.fasta                        
    raw_data = os.path.splitext(sequence)[0] + ".vienna"                        
    VIENNA_PARAMS = os.path.join(os.environ["VIENNA"], "dna_mathews2004.par")

    v_energy = "Undefined"                                                        
    try:                                                                        
        with open(sequence, 'r') as fhin:                                       
            with open(raw_data, 'w') as fhout:
                v_result = subprocess.call(["RNAfold", "--gquad", "--noconv", "--paramFile", VIENNA_PARAMS], stdin=fhin, stdout=fhout)
                ps_file = glob.glob("*.ps")                                     
                if ps_file[0]:                                                  
                    v_result = subprocess.call(["ps2pdf", ps_file[0]])            
    except IOError:                                                             
        pass                                                                    
    else:                                                                       
        with open(raw_data, 'r') as fh:
            for idx, line in enumerate(fh):                                     
                if idx == 2:                                                    
                    v_energy = line.split()[-1].strip("(").strip(")")             
    return v_energy

def fold_seq(input_args):
    """
    list of arguments the first is the filename i.e. 00001.fasta
    and the second is the directory to work in
    """
    os.chdir(input_args[1])
    vienna_energy = run_vienna(input_args[0])
    
    mfold_arg1="SEQ=" + input_args[0]
    mfold_arg2="NA=DNA"
    with open(os.devnull, "w") as fnull:
        #subprocess.check_call(["mfold", mfold_arg1, mfold_arg2], stderr=subprocess.STDOUT, stdout=fnull)
        subprocess.check_call(["mfold", mfold_arg1, mfold_arg2], stderr=fnull, stdout=fnull)
        #p = subprocess.Popen(["mfold", mfold_arg1, mfold_arg2], stdout=fnull, stderr=fnull)
        #p.communicate()

    seq_base_name = os.path.splitext(input_args[0])[0]
    try:
        with open(seq_base_name + ".ct", 'r') as fh:
            for line in fh:
                energy = line.split()[3]
                break
    except IOError:
        energy = "Undefined"            

    targz_name = tarfile.open(input_args[0] + ".tgz",'w:gz')
    for a_file in os.listdir(input_args[1]):
        if a_file.split(".")[-1] != "tgz":
            targz_name.add(a_file)

    targz_name.close()
    for a_file in os.listdir(input_args[1]):
        if a_file.split(".")[-1] != "tgz":
            os.remove(a_file)

    with open("energy.dat", "w") as e_out:
        e_out.write(input_args[0] + " " + energy + " " + vienna_energy)


def driver(scratch_dir, fasta_file_name):

    # Clear existing directories
    try:
        shutil.rmtree(scratch_dir)
    except OSError:
        pass

    os.mkdir(scratch_dir)
    # Create the pool of workers 
    #TODO use celery eventually?
    worker_pool = Pool(NUM_PROCESSORS)
    jobs = []
    # Create workdirs for cpus
    cpu_workdirs = []
    for i in range(NUM_PROCESSORS):
        cpu_workdir = os.path.join(scratch_dir, "cpu_" + str(i))
        os.mkdir(cpu_workdir)
        cpu_workdirs.append(cpu_workdir)

    with open(fasta_file_name, "rU") as handle:
        seq_segments = []
        for record in SeqIO.parse(handle, "fasta"):
            num_seq = len(record.seq)
            print record.id
            print "num seq ", num_seq
            print "divid by window size %f", float(num_seq / SLIDING_WINDOW_SIZE)
            print "module by window size %f", float(num_seq % SLIDING_WINDOW_SIZE)
            seq_segments =  [ record.seq[i:i+FOLD_SIZE] for i in range(0, num_seq, SLIDING_WINDOW_SIZE)]

        print len(seq_segments)


    cpu = 0
    #for i in [0,10001,10003, 100002, 100003, 100004, 10005, 100006]:
    for i in range(len(seq_segments)):
        #print seq_segments[i], len(seq_segments[i])
        prefix = ""
        if i < 10:
            prefix = "0000" + str(i)
        elif i < 100:
            prefix = "000" + str(i)
        elif i < 1000:
            prefix = "00" + str(i)
        elif i < 10000:
            prefix = "0" + str(i)
        else:
            prefix = "" + str(i)

        dir_check = (float(num_seq / SLIDING_WINDOW_SIZE) / NUM_PROCESSORS) * (cpu + 1)

        if i < dir_check:
            pass
        else:
            cpu += 1
        try:
            current_write_dir = cpu_workdirs[cpu]
        except IndexError:
            print "skipping i, cpu", i, cpu
            pass
        else:
            seq_dir = os.path.join(current_write_dir, prefix)
            os.mkdir(seq_dir)
            seq_file = os.path.join(seq_dir, prefix + ".fasta")
            short_record = SeqRecord(seq_segments[i], 'segment_%s' % prefix, '', '')
            output_handle = open(seq_file, "w")
            SeqIO.write(short_record, output_handle, "fasta")
            output_handle.close()
            jobs.append([prefix + ".fasta", seq_dir])

    print "folding", jobs

    worker_pool.map(fold_seq, jobs)

def get_energy_simple(scratch_dir, fasta_file_name):
    logging.basicConfig(filename='build_master_energy.log', level=logging.DEBUG)
    logging.debug(time.strftime("%H:%M:%S"))

    # Clear existing directories
    cpu_workdirs = []
    for i in range(NUM_PROCESSORS):
        cpu_workdir = os.path.join(scratch_dir, "cpu_" + str(i))
        cpu_workdirs.append(cpu_workdir)
    with open(fasta_file_name, "rU") as handle:
        seq_segments = []
        for record in SeqIO.parse(handle, "fasta"):
            num_seq = len(record.seq)
            logging.debug("record id %s"%record.id)
            logging.debug("num seq %s"%num_seq)
            logging.debug("Divide by window size %f"%float(num_seq / SLIDING_WINDOW_SIZE))
            logging.debug("Modulo by window size %f"%float(num_seq % SLIDING_WINDOW_SIZE))
            seq_segments = [record.seq[i:i+FOLD_SIZE] for i in range(0, num_seq, SLIDING_WINDOW_SIZE)]

        logging.debug("Seq segment length : %s" % len(seq_segments))

    cpu = 0
    master_list = []
    for i in range(len(seq_segments)):

        if i < 10:
            prefix = "0000" + str(i)
        elif i < 100:
            prefix = "000" + str(i)
        elif i < 1000:
            prefix = "00" + str(i)
        elif i < 10000:
            prefix = "0" + str(i)
        else:
            prefix = "" + str(i)

        dir_check = (float(num_seq / SLIDING_WINDOW_SIZE) / NUM_PROCESSORS) * (cpu + 1)
        if i < dir_check:
            pass
        else:
            cpu += 1
        try:
            current_write_dir = cpu_workdirs[cpu]
        except IndexError:
            logging.debug("Index error, skipping %s, for cpu %s"%(i, cpu))
            pass
        else:

            seq_dir = os.path.join(current_write_dir, prefix)
            ename = os.path.join(seq_dir, "energy.dat")
            try:
                with open(ename, 'r') as fh:
                    for line in fh: 
                        master_list.append(line)
                        break
            except IOError:
                logging.debug("Could not read energy file,%s" % (ename))
                pass

    master_energy_filename = "master_energy.dat"
    with open(master_energy_filename, 'w') as fh:
        for elem in master_list:
            fh.write(elem)
            fh.write('\n')

def main():
    args = parse_it()
    scratch_dir = os.path.join(os.getcwd(), args.scratch_dir)
    driver(scratch_dir, args.fasta_file_name)
    get_energy_simple(scratch_dir, args.fasta_file_name)

if __name__ == '__main__':
    main()
