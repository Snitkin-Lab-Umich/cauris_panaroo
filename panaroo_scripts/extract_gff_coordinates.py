
import argparse
import gffutils as gff

def extract_gff_coordinates(gff_path, output_path):
    # take a gff file as input, and extract the coordinates of each feature
    # write out a table of feature ID, chromosome, start, end
    gff_db = gff.create_db(gff_path,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    with open(output_path, 'w') as out_f:
        out_f.write('FeatureID\tChromosome\tStart\tEnd\n')
        for feature in gff_db.all_features():
            feature_id = feature.id
            chromosome = feature.seqid
            start = feature.start
            end = feature.end
            out_f.write(f'{feature_id}\t{chromosome}\t{start}\t{end}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide a path to a gff file.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to an output file.''',
        required=True
        )
    args = parser.parse_args()
    extract_gff_coordinates(args.gff, args.output)

if __name__ == '__main__':
    main()
