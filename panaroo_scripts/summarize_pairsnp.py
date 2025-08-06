import os
import argparse
import pandas as pd

# summarizes the results of pairsnp outputs

def get_clades(qc_file):
    cdata = {'Clade I': [], 'Clade III': [], 'Clade IV': []}
    with open(qc_file, 'r') as fh:
        cols = fh.readline().strip().split('\t')
        clade_pos = cols.index('auriclass_clade')
        isolate_pos = cols.index('Sample')
        for line in fh:
            line2 = line.strip().split('\t')
            clade = line2[clade_pos]
            isolate = line2[isolate_pos]
            if clade in cdata:
                cdata[clade].append(isolate)
    return (cdata['Clade I'], cdata['Clade III'], cdata['Clade IV'])


def get_means_normalized_v2(input_file, c1, c3, c4):
    # import the pairsnp output file as a pandas DataFrame and calculate means
    df = pd.read_csv(input_file, sep='\t', header=None, names=['isolate1', 'isolate2', 'distance','gene_length','normalized_distance'])
    df['isolate1'] = df['isolate1'].str.split(';').str[0]
    df['isolate2'] = df['isolate2'].str.split(';').str[0]
    # all comparisons
    all_mean = round(df['normalized_distance'].mean(),5)
    # within clade
    c1_within = df[df['isolate1'].isin(c1) & df['isolate2'].isin(c1)]
    c3_within = df[df['isolate1'].isin(c3) & df['isolate2'].isin(c3)]
    c4_within = df[df['isolate1'].isin(c4) & df['isolate2'].isin(c4)]
    mean_c1 = round(c1_within['normalized_distance'].mean(), 5) if not c1_within.empty else -1
    mean_c3 = round(c3_within['normalized_distance'].mean(), 5) if not c3_within.empty else -1
    mean_c4 = round(c4_within['normalized_distance'].mean(), 5) if not c4_within.empty else -1
    # between clades
    c1_to_c3 = df[(df['isolate1'].isin(c1) & df['isolate2'].isin(c3)) | (df['isolate1'].isin(c3) & df['isolate2'].isin(c1))]
    c1_to_c4 = df[(df['isolate1'].isin(c1) & df['isolate2'].isin(c4)) | (df['isolate1'].isin(c4) & df['isolate2'].isin(c1))]
    c3_to_c4 = df[(df['isolate1'].isin(c3) & df['isolate2'].isin(c4)) | (df['isolate1'].isin(c4) & df['isolate2'].isin(c3))]
    mean_c13 = round(c1_to_c3['normalized_distance'].mean(),5) if not c1_to_c3.empty else -1
    mean_c14 = round(c1_to_c4['normalized_distance'].mean(),5) if not c1_to_c4.empty else -1
    mean_c34 = round(c3_to_c4['normalized_distance'].mean(),5) if not c3_to_c4.empty else -1
    # gene length
    gene_length = df['gene_length'].iloc[0]
    return (all_mean, gene_length, mean_c1, mean_c3, mean_c4, mean_c13, mean_c14, mean_c34)


def get_means_v2(input_file, c1, c3, c4):
    # import the pairsnp output file as a pandas DataFrame and calculate means
    df = pd.read_csv(input_file, sep='\t', header=None, names=['isolate1', 'isolate2', 'distance'])
    df['isolate1'] = df['isolate1'].str.split(';').str[0]
    df['isolate2'] = df['isolate2'].str.split(';').str[0]
    # all comparisons
    all_mean = round(df['distance'].mean(),2)
    # within clade
    c1_within = df[df['isolate1'].isin(c1) & df['isolate2'].isin(c1)]
    c3_within = df[df['isolate1'].isin(c3) & df['isolate2'].isin(c3)]
    c4_within = df[df['isolate1'].isin(c4) & df['isolate2'].isin(c4)]
    mean_c1 = round(c1_within['distance'].mean(), 2) if not c1_within.empty else -1
    mean_c3 = round(c3_within['distance'].mean(), 2) if not c3_within.empty else -1
    mean_c4 = round(c4_within['distance'].mean(), 2) if not c4_within.empty else -1
    # between clades
    c1_to_c3 = df[(df['isolate1'].isin(c1) & df['isolate2'].isin(c3)) | (df['isolate1'].isin(c3) & df['isolate2'].isin(c1))]
    c1_to_c4 = df[(df['isolate1'].isin(c1) & df['isolate2'].isin(c4)) | (df['isolate1'].isin(c4) & df['isolate2'].isin(c1))]
    c3_to_c4 = df[(df['isolate1'].isin(c3) & df['isolate2'].isin(c4)) | (df['isolate1'].isin(c4) & df['isolate2'].isin(c3))]
    mean_c13 = round(c1_to_c3['distance'].mean(),2) if not c1_to_c3.empty else -1
    mean_c14 = round(c1_to_c4['distance'].mean(),2) if not c1_to_c4.empty else -1
    mean_c34 = round(c3_to_c4['distance'].mean(),2) if not c3_to_c4.empty else -1
    return (all_mean, mean_c1, mean_c3, mean_c4, mean_c13, mean_c14, mean_c34)

def get_means(input_file, c1, c3, c4):
    # take a pairsnp output file and return both the overall mean and the means between clades
    all_distances = []
    c1_to_c3 = []
    c1_to_c4 = []
    c3_to_c4 = []
    with open(input_file, 'r') as fh:
        for line in fh:
            line2 = line.strip().split('\t')
            isolate1 = line2[0].split(';')[0]
            isolate2 = line2[1].split(';')[0]
            distance = int(line2[2])
            all_distances.append(distance)
            if (isolate1 in c1 and isolate2 in c3) or (isolate1 in c3 and isolate2 in c1):
                c1_to_c3.append(distance)
            if (isolate1 in c1 and isolate2 in c4) or (isolate1 in c4 and isolate2 in c1):
                c1_to_c4.append(distance)
            if (isolate1 in c3 and isolate2 in c4) or (isolate1 in c4 and isolate2 in c3):
                c3_to_c4.append(distance)
    mean_all = sum(all_distances) / len(all_distances) if all_distances else -1
    mean_c13 = sum(c1_to_c3) / len(c1_to_c3) if c1_to_c3 else -1
    mean_c14 = sum(c1_to_c4) / len(c1_to_c4) if c1_to_c4 else -1
    mean_c34 = sum(c3_to_c4) / len(c3_to_c4) if c3_to_c4 else -1
    return (mean_all, mean_c13, mean_c14, mean_c34)

def summarize_pairsnp(input_dir, output_file, qc_file, flag_normalized):
    c1,c3,c4 = get_clades(qc_file)
    # c1, c3, and c4 are lists of the isolates in each clade
    # b8441 is not included in the QC file
    c1.append('b8441')
    with open(output_file, 'w') as fh_out:
        if flag_normalized:
            fh_out.write('group\tgene_length\tmean_all\tmean_c1\tmean_c3\tmean_c4\tmean_c13\tmean_c14\tmean_c34\n')
        else:
            fh_out.write('group\tmean_all\tmean_c1\tmean_c3\tmean_c4\tmean_c13\tmean_c14\tmean_c34\n')
        for fname in os.listdir(input_dir):
            if fname.endswith('pairsnp'):
                groupname = fname.split('.pairsnp')[0]
                input_file = os.path.join(input_dir, fname)
                if flag_normalized:
                    mean_all,gene_length,mean_c1,mean_c3,mean_c4,mean_c13,mean_c14,mean_c34 = get_means_normalized_v2(input_file, c1, c3, c4)
                    fh_out.write(f'{groupname}\t{gene_length}\t{mean_all}\t{mean_c1}\t{mean_c3}\t{mean_c4}\t{mean_c13}\t{mean_c14}\t{mean_c34}\n')
                else:
                    mean_all,mean_c1,mean_c3,mean_c4,mean_c13,mean_c14,mean_c34 = get_means_v2(input_file, c1, c3, c4)
                    fh_out.write(f'{groupname}\t{mean_all}\t{mean_c1}\t{mean_c3}\t{mean_c4}\t{mean_c13}\t{mean_c14}\t{mean_c34}\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a pairsnp output directory. Files should all end in .pairsnp and be in the -s format''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide an output file for the summary.''',
        required=True
        )
    parser.add_argument(
        '--qc','-q',type=str,
        help='''Provide a QC file to determine which isolates belong to which clades.''',
        required=True
        )
    parser.add_argument(
        '--normalized','-n',action='store_true',
        help='''(Optional) Use this flag if you normalized the pairsnp output files.''',
        )
    args = parser.parse_args()
    summarize_pairsnp(args.input, args.output, args.qc, args.normalized)


if __name__ == '__main__':
    main()
