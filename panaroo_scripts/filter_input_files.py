import os
import argparse
import pandas as pd
import shutil

def filter_fastas(input_dir, filter_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filter_file = pd.read_csv(filter_file, sep=',')
    files_to_keep = filter_file['Sample'].tolist()
    print(f'Keeping {len(files_to_keep)} files based on filter file.')
    # for each file in the input directory, if it contains the sample name, copy it to the output directory
    for filename in os.listdir(input_dir):
        for sample in files_to_keep:
            if sample in filename:
                input_path = os.path.join(input_dir, filename)
                output_path = os.path.join(output_dir, filename)
                shutil.copyfile(input_path, output_path)
                break
    print(f'All files copied to {output_dir}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the directory of protein fastas to be filtered.''',
        default=None
        )
    parser.add_argument(
        '--filter_file','-f',type=str,
        help='''Provide the path to the file that contains the fasta files you want to keep.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output directory.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.input):
        print(f'Could not find input file at {args.input}')
    filter_fastas(args.input, args.filter_file, args.output)

if __name__ == '__main__':
    main()