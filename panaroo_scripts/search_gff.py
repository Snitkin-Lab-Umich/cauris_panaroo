import os
import argparse
import gffutils as gff

def find_cds_locations(input_path, output_path, gff_dir):
    # search through a table containing isolate names and cds names
    # add the location of each cds to the table
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        next(infile)
        _ = outfile.write('isolate\tcds\tlocation\n')
        isolate_dict = {}
        # read through the file and record each cds for each isolate
        for line in infile:
            line = line.strip().split(',')
            isolate = line[0]
            cds = line[1]
            if isolate not in isolate_dict:
                isolate_dict[isolate] = []
            isolate_dict[isolate].append(cds)
        # for each isolate, open the gff file and locate each cds
        for isolate in isolate_dict:
            gff_path = os.path.join(gff_dir, f'{isolate}.gff')
            if not os.path.exists(gff_path):
                print(f'Could not find gff file for {isolate} at {gff_path}')
                quit(1)
            db = gff.create_db(gff_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
            for cds in isolate_dict[isolate]:
                try:
                    feature = db[cds]
                    location = feature.chrom
                    outfile.write(f'{isolate}\t{cds}\t{location}\n')
                except:
                    #print(f'Could not find {cds} in {gff_path}')
                    outfile.write(f'{isolate}\t{cds}\tNA\n')
            

def find_cds_concentration(input_path, output_path, gff_dir):
    # search through a table containing isolate names and rare cds names
    # for each isolate in the list, count the number of rare cds features present on each contig
    # start by generating a dictionary of isolate to list of rare cds features
    with open(input_path, 'r') as infile:
        next(infile)
        isolate_dict = {}
        for line in infile:
            line = line.strip().split(',')
            isolate = line[0]
            cds = line[1]
            if isolate not in isolate_dict:
                isolate_dict[isolate] = []
            isolate_dict[isolate].append(cds)
    # for each isolate, open the gff file and count the number of rare and total cds features on each contig
    scaffold_dict = {}
    for isolate in isolate_dict:
        rare_cds_list = isolate_dict[isolate]
        scaffold_dict[isolate] = {}
        gff_path = os.path.join(gff_dir, f'{isolate}.gff')
        if not os.path.exists(gff_path):
            print(f'Could not find gff file for {isolate} at {gff_path}')
            quit(1)
        # get the total size of every scaffold in the gff file
        with open(gff_path, 'r') as gff_file:
            for line in gff_file:
                if line.startswith('##'):
                    if line.startswith('##sequence-region'):
                        line = line.strip().split()
                        scaffold_name = line[1]
                        scaffold_length = int(line[3]) - int(line[2]) + 1
                        # int(line[2]) should always be 1, but just in case
                        if scaffold_name not in scaffold_dict[isolate]:
                            scaffold_dict[isolate][scaffold_name] = {'total':0, 'rare':0, 'length':scaffold_length}
                        else:
                            print(f'Warning: duplicate scaffold {scaffold_name} in {gff_path}')
        db = gff.create_db(gff_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        for feature in db.all_features():
            if feature.featuretype == 'CDS':
                # all features should be CDS, but just in case
                scaffold = feature.chrom
                scaffold_dict[isolate][scaffold]['total'] += 1
                if feature.id in rare_cds_list:
                    scaffold_dict[isolate][scaffold]['rare'] += 1
    # write the output file
    with open(output_path, 'w') as outfile:
        outfile.write('isolate\tscaffold\ttotal_cds\trare_cds\tscaffold_length\trare_cds_concentration\n')
        for isolate in scaffold_dict:
            for scaffold in scaffold_dict[isolate]:
                total_cds = scaffold_dict[isolate][scaffold]['total']
                rare_cds = scaffold_dict[isolate][scaffold]['rare']
                length = scaffold_dict[isolate][scaffold]['length']
                concentration = rare_cds / total_cds if total_cds > 0 else 0
                outfile.write(f'{isolate}\t{scaffold}\t{total_cds}\t{rare_cds}\t{length}\t{concentration}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the input table.''',
        default=None
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the path to the directory of gff files.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output table.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.input) or not os.path.exists(args.gff):
        print(f'Could not find input file at {args.input} or gff directory at {args.gff}')
    find_cds_concentration(args.input, args.output, args.gff)

if __name__ == '__main__':
    main()