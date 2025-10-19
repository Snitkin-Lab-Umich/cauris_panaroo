import os
import subprocess
import argparse

def create_bwa_alignment(isolate, input_assembly, input_reads_r1, input_reads_r2, output_dir, run_depth=False):
    # create a bwa alignment of the reads to the assembly
    # create a folder in the output directory for the sample
    sample_output_dir = os.path.join(output_dir, 'alignments', isolate)
    if not os.path.exists(sample_output_dir):
        os.makedirs(sample_output_dir)
    # start by copying the assembly to the output directory
    assembly_copy = os.path.join(sample_output_dir, f'{isolate}.scaffolds.fa')
    subprocess.run(['cp', input_assembly, assembly_copy])
    # index the assembly with bwa
    subprocess.run(['bwa', 'index', assembly_copy])
    # align the reads to the assembly
    sam_file = os.path.join(sample_output_dir, f'{isolate}_aligned.sam')
    with open(sam_file, 'w') as sam_fh:
        subprocess.run(['bwa', 'mem', '-t', '8', assembly_copy, input_reads_r1, input_reads_r2], stdout=sam_fh)
    # convert the sam file to a bam file
    bam_file = os.path.join(sample_output_dir, f'{isolate}_aligned.bam')
    subprocess.run(['samtools', 'view', '-b', '-o', bam_file, sam_file])
    # sort the bam file
    sorted_bam_file = os.path.join(sample_output_dir, f'{isolate}_aligned_sorted.bam')
    subprocess.run(['samtools', 'sort', bam_file, '-o', sorted_bam_file])
    # index the sorted bam file
    subprocess.run(['samtools', 'index', sorted_bam_file])
    # remove the intermediate files
    os.remove(sam_file)
    os.remove(bam_file)
    if run_depth:
        # create a depth file
        depth_file = os.path.join(sample_output_dir, f'{isolate}_depth.txt')
        with open(depth_file, 'w') as depth_fh:
            subprocess.run(['samtools', 'depth', '-aa', sorted_bam_file], stdout=depth_fh)
    return(sorted_bam_file)


def calculate_scaffold_coverage(depth_file):
    scaffold_coverage = {}
    with open(depth_file, 'r') as depth_fh:
        for line in depth_fh:
            line = line.strip().split('\t')
            scaffold = line[0]
            position = int(line[1])
            coverage = int(line[2])
            if scaffold not in scaffold_coverage:
                scaffold_coverage[scaffold] = {'total_coverage': 0, 'length': 0}
            scaffold_coverage[scaffold]['total_coverage'] += coverage
            scaffold_coverage[scaffold]['length'] += 1
    return(scaffold_coverage)


def align_and_calculate(input_table, assemblies_dir, raw_reads_dir, raw_reads_dir2, output_dir, perform_alignment=True):
    isolates = set()
    scaffold_coverage_data = {}
    # read in each isolate from the input table
    with open(input_table, 'r') as infile:
        next(infile)
        for line in infile:
            isolate = line.strip().split(',')[0]
            if isolate == 'b8441':
                continue
            isolates.add(isolate)
    print(isolates)
    # for each isolate, create a bwa alignment of the reads to the assembly
    for isolate in isolates:
        if perform_alignment:
            assembly_path = os.path.join(assemblies_dir, f'{isolate}.scaffolds.fa')
            reads_r1_path = os.path.join(raw_reads_dir, f'{isolate}_R1_trim_paired.fastq.gz')
            reads_r2_path = os.path.join(raw_reads_dir, f'{isolate}_R2_trim_paired.fastq.gz')
            if not os.path.exists(reads_r1_path):
                reads_r1_path = os.path.join(raw_reads_dir2, f'{isolate}_R1_trim_paired.fastq.gz')
                reads_r2_path = os.path.join(raw_reads_dir2, f'{isolate}_R2_trim_paired.fastq.gz')
            if not all([os.path.exists(x) for x in [assembly_path, reads_r1_path, reads_r2_path]]):
                print(f'Could not find all input files for {isolate}')
                quit(1)
            sorted_bam = create_bwa_alignment(isolate, assembly_path, reads_r1_path, reads_r2_path, output_dir, run_depth=True)
            print(f'Created alignment for {isolate} at {sorted_bam}')
        # determine the coverage of each scaffold present in the depth file
        depth_file = os.path.join(output_dir, 'alignments', isolate, f'{isolate}_depth.txt')
        scaffold_coverage = calculate_scaffold_coverage(depth_file)
        scaffold_coverage_data[isolate] = scaffold_coverage
    # write the output file
    output_path = os.path.join(output_dir, 'scaffold_coverage.tsv')
    with open(output_path, 'w') as outfile:
        outfile.write('isolate\tscaffold\tscaffold_length\ttotal_coverage\taverage_coverage\n')
        for isolate in scaffold_coverage_data:
            for scaffold in scaffold_coverage_data[isolate]:
                length = scaffold_coverage_data[isolate][scaffold]['length']
                total_cov = scaffold_coverage_data[isolate][scaffold]['total_coverage']
                avg_cov = total_cov / length if length > 0 else 0
                outfile.write(f'{isolate}\t{scaffold}\t{length}\t{total_cov}\t{avg_cov:.2f}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the input table.''',
        default=None
        )
    parser.add_argument(
        '--assemblies','-a',type=str,
        help='''Provide the path to the directory of assemblies.''',
        default=None
        )
    parser.add_argument(
        '--raw-reads','-r',type=str,
        help='''Provide the path to the directory of raw reads.''',
        default=None
        )
    parser.add_argument(
        '--raw-reads-extra','-re',type=str,
        help='''Provide the path to a second directory of raw reads.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output directory.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.input) or not os.path.exists(args.assemblies) or not os.path.exists(args.raw_reads):
        print(f'Could not find input directories')
        quit(1)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    align_and_calculate(args.input, args.assemblies, args.raw_reads, args.raw_reads_extra, args.output, perform_alignment=True)

if __name__ == '__main__':
    main()