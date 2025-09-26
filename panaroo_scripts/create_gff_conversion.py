import os
import argparse
import gffutils as gff

# requires biopython
# take an input GFF and return the old locus tag, new locus tag, and the cds name

def get_gff_data(gff_file):
    all_locus_tags = set()
    lt_cds = {}
    lt_oldlt = {}
    gffdb = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    for feature in gffdb.all_features():
        if feature.featuretype == 'CDS':
            cds_name = feature.attributes.get('ID', ['NA'])[0]
            #locus_tag = feature.attributes.get('locus_tag', ['NA'])[0]
            product = feature.attributes.get('product', ['NA'])[0]
            lt_cds[locus_tag] = [cds_name, product]
        if feature.featuretype == 'gene':
            #locus_tag = feature.attributes.get('locus_tag', ['NA'])[0]
            old_locus_tag = feature.attributes.get('old_locus_tag', ['NA'])[0]
            lt_oldlt[locus_tag] = old_locus_tag
        locus_tag = feature.attributes.get('locus_tag', ['NA'])[0]
        all_locus_tags.add(locus_tag)
    return all_locus_tags, lt_cds, lt_oldlt


def write_gff_data(gff_file, output_file):
    all_locus_tags, lt_cds, lt_oldlt = get_gff_data(gff_file)
    with open(output_file, 'w') as fh_out:
        fh_out.write('locus_tag\tcds_name\tproduct\told_locus_tag\n')
        for lt in sorted(all_locus_tags):
            cds_name, product = lt_cds.get(lt, ['NA', 'NA'])
            old_lt = lt_oldlt.get(lt, 'NA')
            fh_out.write(f'{lt}\t{cds_name}\t{product}\t{old_lt}\n')

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a GFF file.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide an output file for the summary.''',
        required=True
        )
    args = parser.parse_args()
    write_gff_data(args.input, args.output)

if __name__ == '__main__':
    main()
