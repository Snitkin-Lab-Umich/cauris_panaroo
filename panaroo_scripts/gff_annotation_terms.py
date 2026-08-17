import os
import argparse
import gffutils as gff
import pandas as pd

def get_gff_annotation_terms(gff_file, annotation_types, feature_type, output_file):
    # take a gff file as input, as well as a list of annotation types to search for (GO, IPS, PFAM, etc.)
    # return a table with one row for each gene name, with each of the specified terms as an additional column
    gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    annotation_dict = {'GeneName':[], **{annotation_type:[] for annotation_type in annotation_types}} 
    # initialize a dictionary to hold the annotation data, with a key for gene name and a key for each annotation type
    for feature in gff_db.features_of_type(feature_type):
        gene_name = feature.attributes.get('ID', ['NA'])[0]
        annotation_dict['GeneName'].append(gene_name)
        # annotations are stored either in Dbxref (for InterProScan and such) or Note (for Eggnog)
        gene_annotation_list = feature.attributes.get('Dbxref', ['']) + feature.attributes.get('Note', ['']) + feature.attributes.get('note', [''])
        # for each annotation type, search through this list to find any matching annotations, and add them to the dict
        for annotation_type in annotation_types:
            annotation_value = []
            for annotation_entry in gene_annotation_list:
                if annotation_type in annotation_entry:
                    annotation_value.append(annotation_entry)
            if annotation_value == []:
                annotation_dict[annotation_type].append('NA')
            else:
                annotation_dict[annotation_type].append(','.join(annotation_value))
    annotation_df = pd.DataFrame(annotation_dict)
    annotation_df.to_csv(output_file, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the path to the input gff file.''',
        default=None
        )
    parser.add_argument(
        '--feature_type','-ft',choices=['CDS','gene','mRNA'],type=str,
        help='''Provide the feature type to search for.''',
        default='mRNA'
        )
    parser.add_argument(
        '--annotation_types','-at',type=str,nargs='+',
        help='''Provide a list of annotation types to include. Make sure this matches the format of the annotations in the gff file.''',
        default=None
        )
    parser.add_argument(
        '--output_file','-o',type=str,
        help='''Provide the path to the output file.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.gff):
        print(f'Could not find gff file at {args.gff}')
        quit(1)
    get_gff_annotation_terms(args.gff, args.annotation_types, args.feature_type, args.output_file)

if __name__ == '__main__':
    main()