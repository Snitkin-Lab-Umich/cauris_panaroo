import os
import argparse
import gffutils as gff
from Bio import SeqIO

def mask_telomeres(input_fasta, input_gff, bases_masked, output_fasta, feature_buffer=2000):
    # take a gff file and a fasta file
    # mask the first and last bases of each scaffold with N's (hard masking) based on the number of specified bases_masked
    # do not mask nucleotides that are within 2 kb of any annotated feature in the gff file
    out_fasta_records = []
    db = gff.create_db(input_gff, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_length = len(record.seq)
        mask_start = bases_masked
        mask_end = record_length - bases_masked
        for feature in db.region(seqid=record.id, featuretype=('gene','CDS','mRNA','exon')):
            # note that gff files are 1-based
            feature_start = feature.start - feature_buffer - 1
            feature_end = feature.end + feature_buffer - 1
            if feature_start < mask_start:
                print(f'Found {feature.id} at {feature_start}')
                mask_start = feature_start
            if feature_end > mask_end:
                print(f'Found {feature.id} at {feature_end}')
                mask_end = feature_end
        if mask_start < 0:
            mask_start = 0
        if mask_end > record_length:
            mask_end = record_length 
        if mask_end < mask_start:
            print(f'Error: scaffold {record.id} is too short to mask {bases_masked} bases.')
            quit(1)
        print(f'{record.id}: masking from 0 to {mask_start} and from {mask_end} to {record_length}')
        new_seq = 'N' * mask_start + record.seq[mask_start:(mask_end+1)] + 'N' * (record_length - mask_end - 1)
        if len(new_seq) != record_length:
            print(f'Masked sequence length error: {len(new_seq)} != {record_length}')
            quit(1)
        new_record = record
        new_record.seq = new_seq
        out_fasta_records.append(new_record)
    SeqIO.write(out_fasta_records, output_fasta, 'fasta')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the input fasta file.''',
        default=None
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the path to the input gff file.''',
        default=None
        )
    parser.add_argument(
        '--bases_masked','-b',type=int,
        help='''Provide the number of bases to mask at the beginning and end of each scaffold.''',
        default=5000
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output fasta file.''',
        default=None
        )
    parser.add_argument(
        '--feature_buffer','-fb',type=int,
        help='''Provide the number of bases to buffer around each feature.''',
        default=2000
        )
    args = parser.parse_args()
    mask_telomeres(args.input, args.gff, args.bases_masked, args.output, args.feature_buffer)

if __name__ == '__main__':
    main()