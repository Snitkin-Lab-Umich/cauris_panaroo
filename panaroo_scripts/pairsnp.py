import os
import subprocess
import argparse

# runs pairsnp on all .aln.fas files in the input directory

def run_pairsnp(input_dir, output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            input_file = os.path.join(input_dir, fname)
            output_file = os.path.join(output_dir, groupname + '.pairsnp')
            with open(output_file, 'w') as fh_out:
                _ = subprocess.call(['pairsnp', '-s', input_file], stdout=fh_out)

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
    args = parser.parse_args()
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    run_pairsnp(args.input, args.output)


if __name__ == '__main__':
    main()
