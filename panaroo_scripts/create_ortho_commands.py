import os
import argparse

def copy_diamond_commands(input_path, output_path):
    linecount = 0
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('diamond blastp'):
                outfile.write(line)
                linecount += 1
    print(f'Copied {linecount} diamond commands from {input_path} to {output_path}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the output from OrthoFinder's -op option.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the new text file that will contain the diamond commands.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.input):
        print(f'Could not find input file at {args.input}')
    copy_diamond_commands(args.input, args.output)

if __name__ == '__main__':
    main()