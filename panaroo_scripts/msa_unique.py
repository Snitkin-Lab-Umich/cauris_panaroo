import os
import argparse
import pandas as pd
from Bio import SeqIO

# runs pairsnp on all .aln.fas files in the input directory

# convert an MSA file to a list of the unique SNPs present in each sequence
def msa_unique(input_file):
    ref = None
    diff_dict = {}
    isolatelist = []
    with open(input_file, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            isolatename = record.id.split(';')[0]
            isolatelist.append(isolatename)
            if ref is None:
                ref = record.seq
                seqlength = len(ref)
                continue
            # compare each sequence to the reference
            diff_dict = get_diffs(ref, record.seq, diff_dict, isolatename)
    return diff_dict, isolatelist, seqlength


def get_diffs(ref, query, diff_dict, isolatename):
    '''Return a dictionary of the differences between two sequences from a MSA. Both must be the same length.'''
    # note that insertions will have the full sequence appended to the dictionary key
    ref = ref + 'E'  # add a dummy character to the end to avoid index errors
    query = query + 'E'
    ins_start, del_start = None, None
    for pos, (char1, char2) in enumerate(zip(ref, query + 'E')):
        diff = None
        if char1 != char2:
            if char1 == '-':
                if ins_start is None:
                    ins_start = pos
                    ins_seq = ''
                ins_seq += char2
                if ref[pos+1] != '-':
                    ins_end = pos
                    diff = f'INS{ins_start+1}-{ins_end+1}_{ins_seq}'
                    ins_start = None
            elif char2 == '-':
                if del_start is None:
                    del_start = pos
                if query[pos+1] != '-':
                    del_end = pos
                    diff = f'DEL{del_start+1}-{del_end+1}'
                    del_start = None
            else:
                diff = f'{char1}{pos+1}{char2}'
        if diff is not None:
            if diff not in diff_dict:
                diff_dict[diff] = [isolatename]
            else:
                diff_dict[diff].append(isolatename)
    return(diff_dict)

def make_df(diff_dict, isolatelist, seqlength, output_file):
    # make a dataframe, where the rows are the isolates and the columns are keys in diff_dict
    data_dict = {'isolate': isolatelist, 'seqlength': [seqlength] * len(isolatelist)}
    for snp in diff_dict:
        data_dict[snp] = [1 if isolate in diff_dict[snp] else 0 for isolate in isolatelist]
    df = pd.DataFrame(data_dict)
    # fix the insert names by removing the full sequence appended to the end
    for name in df.columns:
        if 'INS' in name:
            new_name = name_fixer(name, df.columns)
            df = df.rename(columns={name: new_name})
    df.to_csv(output_file, index=False)

def name_fixer(name, namelist):
    '''Remove the sequence from the end of a variant name and replace it with a unique number.'''
    val = 1
    new_name = f'{name.split("_")[0]}_{val}'
    while new_name in namelist:
        val += 1
        new_name = f'{name.split("_")[0]}_{val}'
    return new_name


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
        help='''Provide an output directory for the unique variant matrices.''',
        required=True
        )
    args = parser.parse_args()
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    for file in os.listdir(args.input):
        if file.endswith('.aln.fas'):
            input_file = os.path.join(args.input, file)
            output_file = os.path.join(args.output, file.replace('.aln.fas', '_var_matrix.csv'))
            diff_dict, isolatelist, seqlength = msa_unique(input_file)
            make_df(diff_dict, isolatelist, seqlength, output_file)

if __name__ == '__main__':
    main()
