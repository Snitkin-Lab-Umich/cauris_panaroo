import os
import argparse
import pandas as pd
from collections import Counter

# summarizes the results of pairsnp outputs


def get_clades(qc_file):
    cdata = {'b8441': 'Clade I'}
    m = pd.read_csv(qc_file, sep='\t')
    for i in m.index:
        if m.loc[i, 'QC_EVALUATION'] == 'PASS':
            clade = m.loc[i, 'auriclass_clade']
            isolate = m.loc[i, 'Sample']
            cdata[isolate] = clade
    return (cdata)

def count_unique_variants(input_file, cdata):
    # import presence/absence matrix as a pandas DataFrame
    df = pd.read_csv(input_file, sep=',', header=0)
    # classify each variant depending on which clades it appears in
    vdict = {'unclassified':0, 'intermediate_all':0, 'minor_clade_variant':0}
    for c in ['c1', 'c3', 'c4']:
        vdict[f'{c}_fixed_specific'] = 0
        vdict[f'{c}_fixed_absent'] = 0
        vdict[f'{c}_intermediate_specific'] = 0
        vdict[f'{c}_intermediate_absent'] = 0
    # paralogs appear to have names such as '_R_SRR19738667', and will thus not appear in cdata
    plogs = [x not in cdata and x.startswith('_') for x in df['isolate']]
    # add paralogs to vdict, then remove them from the dataframe
    vdict['paralogs'] = sum(plogs)
    df = df[[not x for x in plogs]]
    # this makes a separate entry for paralogous genes, treating them the same as genes from separate isolates
    isolate_vector = df['isolate']
    clade_vector = [cdata[i] for i in isolate_vector]
    clade_counts = Counter(clade_vector)
    for variant in df.columns:
        if variant == 'isolate':
            continue
        if variant == 'seqlength':
            seqlength = df[variant].iloc[0]
            continue
        variant_vector = df[variant]
        variant_clades = Counter([clade for clade, variant in zip(clade_vector, variant_vector) if variant == 1])
        # this is a dictionary where the clade is the key and the number of times this isolate appears in the clade is the value
        vtype = variant_classifier(clade_counts, variant_clades)
        vdict[vtype] += 1
    return vdict, seqlength

def variant_classifier(clade_counts, variant_clades, upper_limit=0.95, lower_limit=0.02):
    # classify each variant into one of the following categories:
    # clade-specific presence (variant appears in > 95% of a single clade and < 5% in the others)
    # clade-specific absence (variant appears in > 95% of a two clades and < 5% in the remaining)
    # within-clade variant (variant appears between 5% and 95% in one, two, or three clades)
    # NOTE: variants that appear outside of the three main clades are ignored
    c1_presence = variant_clades.get('Clade I', 0) / clade_counts.get('Clade I', 1)
    c3_presence = variant_clades.get('Clade III', 0) / clade_counts.get('Clade III', 1)
    c4_presence = variant_clades.get('Clade IV', 0) / clade_counts.get('Clade IV', 1)
    cvec = [c1_presence, c3_presence, c4_presence]
    namevec = ['c1', 'c3', 'c4']
    # if exactly one clade has this variant:
    if sum([x > upper_limit for x in cvec]) == 1 and sum([x < lower_limit for x in cvec]) == 2:
        c = [c for c,v in zip(namevec, cvec) if v > upper_limit][0]
        return f'{c}_fixed_specific'
    # if exactly two clades have this variant:
    elif sum([x > upper_limit for x in cvec]) == 2 and sum([x < lower_limit for x in cvec]) == 1:
        c = [c for c,v in zip(namevec, cvec) if v < lower_limit][0]
        return f'{c}_fixed_absent'
    # if the variant has intermediate presence in exactly one clade:
    elif sum([lower_limit <= x <= upper_limit for x in cvec]) == 1 and sum([x < lower_limit for x in cvec]) == 2:
        c = [c for c,v in zip(namevec, cvec) if lower_limit <= v <= upper_limit][0]
        return f'{c}_intermediate_specific'
    # if the variant has intermediate presence in exactly two clades:
    elif sum([lower_limit <= x <= upper_limit for x in cvec]) == 2 and sum([x < lower_limit for x in cvec]) == 1:
        c = [c for c,v in zip(namevec, cvec) if v < lower_limit][0]
        return f'{c}_intermediate_absent'
    # if the variant has intermediate presence in all three clades:
    elif all([x == 0 for x in cvec]):
        # several variants exist between the smaller clades (II, V, and VI) that I don't want to examine currently
        return 'minor_clade_variant'
    elif sum([lower_limit <= x <= upper_limit for x in cvec]) == 3:
        return 'intermediate_all'
    else:
        return 'unclassified'

def summarize_unique(input_dir, output_file, qc_file):
    header_written = False
    cdata = get_clades(qc_file)
    with open(output_file, 'w') as fh_out:
        for fname in os.listdir(input_dir):
            if fname.endswith('var_matrix.csv'):
                input_file = os.path.join(input_dir, fname)
                vdict, seqlength = count_unique_variants(input_file, cdata)
                if not header_written:
                    header_data = ['group','gene_length'] + sorted(list(vdict.keys()))
                    fh_out.write('\t'.join(header_data) + '\n')
                    header_written = True
                groupname = fname.split('.var_matrix.csv')[0]
                vdict['group'] = groupname
                vdict['gene_length'] = seqlength
                fh_out.write('\t'.join([str(vdict[x]) for x in header_data]) + '\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to an input directory. Files should all end in .var_matrix.csv''',
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
    args = parser.parse_args()
    summarize_unique(args.input, args.output, args.qc)


if __name__ == '__main__':
    main()
