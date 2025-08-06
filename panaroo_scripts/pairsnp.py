import os
import subprocess
import argparse
from Bio import SeqIO

# runs pairsnp on all .aln.fas files in the input directory

def run_pairsnp(input_dir, output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            input_file = os.path.join(input_dir, fname)
            output_file = os.path.join(output_dir, groupname + '.pairsnp')
            with open(output_file, 'w') as fh_out:
                _ = subprocess.call(['pairsnp', '-s', input_file], stdout=fh_out)

def add_gene_lengths(input_dir, output_dir, msa_dir):
    '''Take a directory of pairsnp output files and add the gene lengths to the output files.'''
    for fname in os.listdir(input_dir):
        if fname.endswith('.pairsnp'):
            input_file = os.path.join(input_dir, fname)
            output_file = os.path.join(output_dir, fname.replace('.pairsnp', '.norm_pairsnp'))
            msa_file = os.path.join(msa_dir, fname.replace('.pairsnp', '.aln.fas'))
            with open(msa_file, 'r') as msa_fh:
                for record in SeqIO.parse(msa_fh, 'fasta'):
                    gene_length = len(record.seq)
                    # all genes in the MSA should have the same length
                    break
            with open(input_file, 'r') as fh_in, open(output_file, 'w') as fh_out:
                for line in fh_in:
                    line = line.strip().split('\t')
                    line2 = line + [str(gene_length), str(int(line[-1])/gene_length)]
                    fh_out.write('\t'.join(line2) + '\n')

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
        help='''Provide an output directory for the distance matrices.''',
        required=True
        )
    parser.add_argument(
        '--normalized','-n',type=str,
        help='''(Optional) Provide another output directory for the distance matrices, normalized by gene length.''',
        default=None
        )
    parser.add_argument(
        '--skip_pairsnp','-s',action='store_true',
        help='''(Optional) Skip running pairsnp and assume that the relevant files are already in the output directory.''',
        )
    args = parser.parse_args()
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    if not args.skip_pairsnp:
        run_pairsnp(args.input, args.output)
    if args.normalized is not None:
        if not os.path.isdir(args.normalized):
            os.makedirs(args.normalized)
        add_gene_lengths(args.output, args.normalized, args.input)


if __name__ == '__main__':
    main()
