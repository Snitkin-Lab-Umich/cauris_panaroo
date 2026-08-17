import os
import argparse
import gffutils as gff
from Bio import SeqIO

def extract_single_annotation(feature_name, fasta_file, gff_file):
    # take a gff file and a fasta file, and extract the sequence of a single feature from the fasta file based on the coordinates in the gff file
    # write the extracted sequence to a new fasta file with the feature name as the header
    gff_db = gff.create_db(gff_file, dbfn=":memory:", force=True, keep_order=True, merge_strategy="create_unique", sort_attribute_values=True)
    feature = gff_db[feature_name]
    scaffold = feature.seqid
    start = feature.start
    end = feature.end
    # extract the sequence from the fasta file
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == scaffold:
            seq = record.seq[start-1:end]
            return(seq)
    print(f'Unable to locate {feature_name} in {fasta_file}')
    return(None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--feature','-f',type=str,
        help='''Provide the feature name.''',
        default=None
        )
    parser.add_argument(
        '--fasta','-ff',type=str,
        help='''Provide the fasta file.''',
        default=None
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the gff file.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''(Optional)Provide the output file.''',
        default=None
        )
    args = parser.parse_args()
    o = extract_single_annotation(args.feature, args.fasta, args.gff)
    if o is None:
        print(f'Unable to locate {args.feature} in {args.fasta}')
    else:
        if args.output is None:
            print(f'>{args.feature}\n{o}')
        else:
            with open(args.output, 'w') as fh:
                fh.write(f'>{args.feature}\n{o}\n')


if __name__ == '__main__':
    main()