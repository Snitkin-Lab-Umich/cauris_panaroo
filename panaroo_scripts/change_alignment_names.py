import os
import argparse
from Bio import SeqIO


def change_names_paralogskip(input_dir,output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            with open(input_dir + fname,'r') as fh_in, open(output_dir + groupname + '.fasta','w') as fh_out:
                seen,towrite = set(),[]
                for record in SeqIO.parse(fh_in,'fasta'):
                    isolatename = record.id.split(';')[0]
                    # THIS WILL SKIP ANY PARALOGS
                    if isolatename not in seen:
                        record.id = isolatename
                        record.name = ''
                        record.description = ''
                        towrite.append(record)
                        seen.add(isolatename)
                SeqIO.write(towrite,fh_out,'fasta')

def change_names(input_dir,output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            with open(input_dir + fname,'r') as fh_in, open(output_dir + groupname + '.fasta','w') as fh_out:
                for line in fh_in:
                    if line.startswith('>'):
                        line = line.split(';')[0] + '\n'
                    _ = fh_out.write(line)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a panaroo-msa output directory. Files should all end in .aln.fas''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide an output directory for the renamed fasta files.''',
        required=True
        )
    args = parser.parse_args()
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    change_names_paralogskip(args.input, args.output)


if __name__ == '__main__':
    main()
    ## this script is meant to rename the output files from panaroo-msa, removing the numbers after the semicolon
    ## currently, this will also remove any paralogs from the msa file
    ## obviously, this is not ideal, but it allows it to work with IQTREE
