from Bio import SeqIO
import argparse

def make_subfasta(input_fasta, scaffold, start, end, output_fasta):
    # take an input fasta file and a set of coordinates
    # write a new fasta file containing only the specified region
    # if only a scaffold is provided, write the whole thing
    with open(output_fasta, 'w') as out_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id == scaffold:
                if start is not None and end is not None:
                    # if the end is longer than the length of the record, adjust it to the length of the record
                    if end > len(record.seq):
                        print('Adjusting end position to length of record')
                        end = len(record.seq)
                    # if the start is larger than the length of the record, return an error
                    if start > len(record.seq):
                        print('Start position is larger than length of record! Cannot extract sequence.')
                        quit(1)
                    sub_seq = record.seq[start-1:end]  # fasta is 1-based, python is 0-based
                    sub_record = SeqIO.SeqRecord(sub_seq, id=f'{record.id}_{start}_{end}', description='')
                    SeqIO.write(sub_record, out_fasta, 'fasta')
                else:
                    SeqIO.write(record, out_fasta, 'fasta')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide an input fasta file.''',
        default=None
        )
    parser.add_argument(
        '--scaffold','-s',type=str,
        help='''Provide the name of a the scaffold to extract sequences for.''',
        default=None
        )
    parser.add_argument(
        '--start','-S',type=int,
        help='''Provide the start position of the region to extract sequences for.''',
        default=None
        )
    parser.add_argument(
        '--end','-E',type=int,
        help='''Provide the end position of the region to extract sequences for.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to an output file.''',
        default=None
        )
    args = parser.parse_args()
    make_subfasta(args.input, args.scaffold, args.start, args.end, args.output)


if __name__ == '__main__':
    main()