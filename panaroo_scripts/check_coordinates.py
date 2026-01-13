
from Bio import SeqIO
import argparse

def get_assembly_sequence(assembly_path, chromosome, start, end):
    for record in SeqIO.parse(assembly_path, 'fasta'):
        if record.id == chromosome:
            sequence = record.seq[start-1:end]  # converting to 0-based indexing
            print(f'>{chromosome}:{start}-{end}\n{sequence}')
            return None
    print(f'Chromosome {chromosome} not found in assembly {assembly_path}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--assembly','-a',type=str,
        help='''Provide a path to an assembly.''',
        required=True
        )
    parser.add_argument(
        '--chromosome','-c',type=str,
        help='''Provide the name of the chromosome.''',
        required=True
        )
    parser.add_argument(
        '--start','-s',type=int,
        help='''Provide the start position.''',
        required=True
        )
    parser.add_argument(
        '--end','-e',type=int,
        help='''Provide the end position.''',
        required=True
        )
    args = parser.parse_args()
    get_assembly_sequence(args.assembly, args.chromosome, args.start, args.end)

if __name__ == '__main__':
    main()
