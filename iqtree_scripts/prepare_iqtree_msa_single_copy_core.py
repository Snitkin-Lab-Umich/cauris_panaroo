import os
import argparse
from Bio import SeqIO
import gffutils as gff
import subprocess
import shutil
import pandas as pd

def check_single_copy_core(input_fasta, isolate_list, outgroup_list=[], core_thresh = 0.95, retain_refound = True):
    # take an msa fasta file and check each name
    # make sure each entry consists only of the isolate name from name_list
    # additionally, if any duplicates entries for the same isolate are found, return False
    # also return False if less than [core_thresh] of the isolates from isolate_list are present in the file
    isolates_seen = set()
    name_list = isolate_list + outgroup_list
    for record in SeqIO.parse(input_fasta, "fasta"):
        isolate_name = record.id.split(';')[0]
        if isolate_name.startswith('_R_'):
            if not retain_refound:
                continue
            isolate_name = isolate_name.replace('_R_','')
        if isolate_name not in name_list:
            print(f'Warning: unexpected isolate name {isolate_name} found in {input_fasta}')
        if isolate_name in isolates_seen:
            return False
        isolates_seen.add(isolate_name)
    # check if the number of isolates seen meets the core threshold from name_list
    if len(isolates_seen) < (core_thresh * len(name_list)):
        return False
    else:
        return True

def fix_alignment_names(input_fasta, output_fasta, isolate_list, outgroup_list=[], retain_refound = True):
    # take an msa fasta file and check each name
    # make sure each entry consists only of the isolate name from name_list
    # additionally, if any duplicates entries for the same isolate are found, raise an error
    isolates_seen = set()
    name_list = isolate_list + outgroup_list
    with open(output_fasta, 'w') as output_fh:
        for record in SeqIO.parse(input_fasta, "fasta"):
            isolate_name = record.id.split(';')[0]
            if isolate_name.startswith('_R_'):
                if not retain_refound:
                    continue
                isolate_name = isolate_name.replace('_R_','')
            if isolate_name not in name_list:
                print(f'Warning: unexpected isolate name {isolate_name} found in {input_fasta}')
            if isolate_name in isolates_seen:
                raise ValueError(f'Duplicate isolate name {isolate_name} found in {input_fasta}')
            isolates_seen.add(isolate_name)
            record.id = isolate_name
            record.description = ''
            SeqIO.write(record, output_fh, "fasta")

def get_name_list(qc_file):
    # read in the master qc file and get the sample column
    # return a list
    df = pd.read_csv(qc_file, sep=',')
    return df['Sample'].tolist() + ['b8441']

def add_outgroup_to_alignment(input_fasta_dir, output_fasta_dir, outgroup_fasta, outgroup_gff, ortholog_file, temp_dir, outgroup_name):
    # use an input table to determine which outgroup CDS sequences to add to each msa file in the input directory
    # each msa file has a name that corresponds to a group in the ortholog file
    # use MAFFT to add sequences

    # read in the ortholog file
    ortho_dict = {}
    with open(ortholog_file, 'r') as fh:
        next(fh)
        for line in fh:
            group, ortho = line.strip().split('\t')
            ortho_dict[group] = ortho
    
    # prepare gff and fasta from the outgroup
    outgroup_gff_db = gff.create_db(outgroup_gff, dbfn=':memory:', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
    outgroup_fasta = {record.id:record for record in SeqIO.parse(outgroup_fasta,'fasta')}
    
    # iterate through each msa file in the input directory, making sure the same cds is never used twice
    used_outgroup_cds = set()
    for msa_file in os.listdir(input_fasta_dir):
        if not msa_file.endswith('.aln.fas'):
            continue
        group_name = msa_file.replace('.aln.fas','')
        if group_name not in ortho_dict:
            continue
        outgroup_ortho = f'CGD:{ortho_dict[group_name]}'
        # outgroup_ortho is the CGD ID of the CDS, which is present in the CGD Dbxref attribute
        outgroup_cds = None
        for feature in outgroup_gff_db.features_of_type('CDS'):
            dbxref = feature.attributes.get('Dbxref',[])
            if outgroup_ortho in dbxref:
                # add to the used cds set
                if feature.id in used_outgroup_cds:
                    print(f'Error: outgroup CDS {feature.id} already used in another alignment')
                    quit(1)
                used_outgroup_cds.add(feature.id)
                # get the corresponding sequence from the fasta
                outgroup_cds = outgroup_fasta[feature.seqid][feature.start-1:feature.end]
                break
        if outgroup_cds is None:
            # copy the msa file without any changes
            shutil.copyfile(os.path.join(input_fasta_dir, msa_file), os.path.join(output_fasta_dir, msa_file))
        else:
            # add the outgroup sequence to the alignment using mafft --add
            # first, write only the outgroup sequence to a temporary fasta file
            temp_outgroup_fasta = os.path.join(temp_dir, 'outgroup_temp.fasta')
            outgroup_cds.id = outgroup_name
            outgroup_cds.description = ''
            with open(temp_outgroup_fasta, 'w') as temp_fh:
                SeqIO.write(outgroup_cds, temp_fh, 'fasta')
            # now run mafft --add
            input_msa = os.path.join(input_fasta_dir, msa_file)
            output_msa = os.path.join(output_fasta_dir, msa_file)
            mafft_command = f'mafft --add {temp_outgroup_fasta} --reorder {input_msa} > {output_msa}'
            subprocess.run(mafft_command, shell=True, check=True)

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
        help='''Provide a path to a new directory where the modified alignment files will go.''',
        required=True
        )
    parser.add_argument(
        '--qc_file','-qc',type=str,
        help='''Provide a path to a qc file to get sample names from.''',
        required=True
        )
    parser.add_argument(
        '--outgroup_gff','-ogff',type=str,
        help='''Provide a path to a gff file for the outgroup you want to add.''',
        )
    parser.add_argument(
        '--outgroup_fasta','-ofasta',type=str,
        help='''Provide a path to a fasta file for the outgroup you want to add.''',
        )
    parser.add_argument(
        '--ortholog','-ortho',type=str,
        help='''Provide a path to a two-column file with the group name in the first column and the outgroup's ortholog in the second column.''',
        )
    parser.add_argument(
        '--outgroup_name','-oname',type=str,
        help='''Provide the name of the outgroup.''',
        )
    parser.add_argument(
        '--outgroup_list','-olist',type=str, nargs='*', default=[],
        help='''Provide a list of outgroup names to check the MSA files for. Files that only contain outgroup sequences will be removed.''',
        )
    args = parser.parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if args.ortholog:
        # make a temp dir
        temp_dir = f'temp/iqtree_{args.outgroup_name}_msa_{str(os.getpid())}'
        temp_dir2 = f'temp/iqtree_{args.outgroup_name}_msa_{str(os.getpid())}/msa_fasta'
        os.makedirs(temp_dir2, exist_ok=True)
        add_outgroup_to_alignment(args.input, temp_dir, args.outgroup_fasta, args.outgroup_gff, args.ortholog, temp_dir2, args.outgroup_name)
        args.input = temp_dir
    name_list = get_name_list(args.qc_file)
    for msa_file in os.listdir(args.input):
        if not msa_file.endswith('.aln.fas'):
            continue
        input_msa = os.path.join(args.input, msa_file)
        output_msa = os.path.join(args.output, msa_file)
        if not check_single_copy_core(input_msa, isolate_list = name_list, outgroup_list = args.outgroup_list):
            print(f'Skipping {input_msa}, as it is not a single-copy core gene family')
            continue
        print(f'Processing {input_msa} to {output_msa}')
        fix_alignment_names(input_msa,output_msa,isolate_list = name_list, outgroup_list = args.outgroup_list)


if __name__ == '__main__':
    # input = 'test/original_msa/group_1000.aln.fas'
    # output = 'test/iqtree_msa/group_1000.aln.fas'
    # qc_file = '092525_pangenome_results/subset700_data/master_qc_file_PassOnly_v2.csv'
    # name_list = get_name_list(qc_file)
    # fix_alignment_names(input, output, name_list)
    main()