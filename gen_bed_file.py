#!/usr/bin/python
import os,sys
import argparse
header="""browser position chr21:42300001-51000000
track type=wiggle_0 name="{TYPE}" description="Replication Data" visibility=full autoScale=off viewLimits={MIN}:0.0 color=255,0,0 yLineMark=11.76 yLineOnOff=on priority=10
variableStep chrom=chr21 span=150
"""

CHOICE={1:"MFOLD", 2:"VIENNA"}


def parse_it():
    parser = argparse.ArgumentParser(description='Process a energy.dat file and generate a bed file.')
    parser.add_argument('energy_file',help='The name of the energy file i.e. energy.dat.')
    parser.add_argument('bed_file', help='The name of the output file i.e. abl.bed')
    parser.add_argument('column', type=int, choices=range(1,3), help=("1=mfold 2=vienna"))
    parser.add_argument('start_pos', type=int, help=("The start position of the folded region i.e. 100323"), default=0)
    parser.add_argument('chromosome', help=("The chromosome for the folded region i.e. chr1"))
    
    args = parser.parse_args()
    return args

def gen_bed(args):
    """
    """
    data_selection = args.column
    chrom = args.chromosome
    start_pos = args.start_pos

    min_val = 0.0
    output=[]
    output_simple=[]
    with open(args.energy_file, 'r') as fin:
        for idx,line in enumerate(fin, start=1):
            val = line.strip().split()[data_selection]
            count = "%s"%(idx * 150)
            abs_count = ((idx - 1) * 150)
            absolute_pos = abs_count + start_pos
            absolute_end_pos = abs_count + start_pos + 300
            if val == "Undefined":
                output.append("%s %s\n"%(count, 0))
                output_simple.append("%s\t%s\t%s\t%s\n"%(chrom,
                                                        absolute_pos, 
                                                        absolute_end_pos, val))
            else:
                output.append("%s %s\n"%(count, val))
                output_simple.append("%s\t%s\t%s\t%s\n"%(chrom,
                                                        absolute_pos, 
                                                        absolute_end_pos, val))
                if float(val) < min_val:
                    min_val = float(val)

    with open(args.bed_file, 'w') as fout:
        fout.write(header.format(TYPE=CHOICE[data_selection], MIN=min_val))
        for d in output:
            fout.write(d)

    with open("simple_%s" % args.bed_file, 'w') as fout:
        for d in output_simple:
            fout.write(d)

    print "min --> ", min_val

def main():
    args = parse_it()
    gen_bed(args)

if __name__ == '__main__':
    main()
