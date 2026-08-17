import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
import gffutils as gff

# the goal of this script is to take a list of gene families from a pangenome and make the files needed to visualize presence and absence via gggenomes
# for each gene family, find three isolates with and three without the gene family
# use the neighbor database to find the correct scaffold for each of these isolates
# subset the fasta and gff files to include only these scaffolds
# align all of the scaffolds using minimap2 and output a paf file


def select_isolates(gene_family, pan, filter_list, good_isolate_prefixes = ['b8441', 'Chi_', 'UM_']):
    gf_series = pan[pan['Gene'] == gene_family].squeeze()
    # change any values that contain the string 'refound' or 'pseudo' to empty strings
    gf_series = gf_series.apply(lambda x: '' if 'refound' in x or 'pseudo' in x else x)
    isolates_with_gf = gf_series[gf_series != ''].index.tolist()
    isolates_without_gf = gf_series[gf_series == ''].index.tolist()
    # subset both lists to include only isolates in the filter list
    isolates_with_gf = [isolate for isolate in isolates_with_gf if isolate in filter_list]
    isolates_without_gf = [isolate for isolate in isolates_without_gf if isolate in filter_list]
    # sort both lists to start with good isolates
    isolates_with_gf = sorted(isolates_with_gf, key = lambda x: 0 if any([x.startswith(prefix) for prefix in good_isolate_prefixes]) else 1)
    isolates_without_gf = sorted(isolates_without_gf, key = lambda x: 0 if any([x.startswith(prefix) for prefix in good_isolate_prefixes]) else 1)
    return isolates_with_gf[:min(3, len(isolates_with_gf))], isolates_without_gf[:min(3, len(isolates_without_gf))]

def find_scaffold(gene_family, neighbor_df):
    # subset to the gene family
    n_series = neighbor_df[neighbor_df['GeneFamily'] == gene_family].squeeze()
    scaffold = n_series['Scaffold']
    if scaffold == 'absent':
        scaffold = n_series['Neighbor1_Scaffold']
    return scaffold

def find_position(gene_family, neighbor_df):
    n_series = neighbor_df[neighbor_df['GeneFamily'] == gene_family].squeeze()
    gene_name = n_series['GeneName']
    scaffold = n_series['Scaffold']
    start = n_series['Start']
    end = n_series['End']
    target_data = [n_series['GeneName'], n_series['Scaffold'], n_series['Start'], n_series['End'], 'TARGET']
    neighbor1_data = [n_series['Neighbor1_GeneName'], n_series['Neighbor1_Scaffold'], n_series['Neighbor1_Start'], n_series['Neighbor1_End'], 'NEIGHBOR']
    neighbor2_data = [n_series['Neighbor2_GeneName'], n_series['Neighbor2_Scaffold'], n_series['Neighbor2_Start'], n_series['Neighbor2_End'], 'NEIGHBOR']
    return [target_data, neighbor1_data, neighbor2_data]

def trim_cds_name(cds_name):
    if cds_name.startswith('cds-'):
        return cds_name.replace('cds-', '', 1)
    elif cds_name.endswith('-T1.cds'):
        return cds_name.replace('-T1.cds', '')
    else:
        return cds_name


def subset_fasta_and_gff_v2(original_fasta_file, original_gff_file, scaffold, start, end, new_fasta_file, new_gff_file, target_gene_name, neighbor_gene_names, constant_gene_name,margin = 5000):
    target_gene_id = trim_cds_name(target_gene_name)
    neighbor_gene_ids = [trim_cds_name(name) for name in neighbor_gene_names]
    with open(new_fasta_file, 'w') as fasta_out, open(new_gff_file, 'w') as gff_out:
        # first, determine if the target gene is on the forward or reverse strand
        gff_db = gff.create_db(original_gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
        # the neighbor gene names should exactly match the ID attribute of a feature in the gff
        # for now, assume that neighbor 1 is usually the same gene
        # invert the coordinates if neighbor 1 is on the reverse strand
        if constant_gene_name is not None:
            target_strand = gff_db[constant_gene_name].strand
        else:
            target_strand = '+'
        # modify the fasta file to contain only the specified scaffold, start, and stop (with margin)
        # also reverse-complement the sequence if the target gene is on the reverse strand
        for record in SeqIO.parse(original_fasta_file, 'fasta'):
            if record.id == scaffold:
                # record the original length of the scaffold for gff coordinate updating
                original_length = len(record.seq)
                # adjust the start and end to include the margin (if possible)
                start = max(1,start - margin)
                end = min(original_length,end + margin)
                # subset the sequence to the start and end
                subset_seq = record.seq[start:end]
                # reverse-complement the sequence if the target gene is on the reverse strand
                if target_strand == '-':
                    subset_seq = subset_seq.reverse_complement()
                SeqIO.write(SeqIO.SeqRecord(subset_seq, id=record.id, description=''), fasta_out, 'fasta')
        # next, modify the gff file to contain only the features in the specified region
        # updated coordinates as needed, including reverse-complementing if the strand is -
        seq_offset = start - 1
        for feature in gff_db.region(seqid=scaffold, start=start, end=end):
            if feature.end < start or feature.start > end:
                continue
            new_start = feature.start - seq_offset
            new_end = feature.end - seq_offset
            new_start = max(1,new_start)
            new_end = min(end - seq_offset,new_end)
            if new_start >= new_end or new_end <= 1:
                continue
            # invert the coordinates of this feature if the strand is -
            if target_strand == '-':
                seq_length = end - start + 1
                new_start, new_end = seq_length - new_end + 1, seq_length - new_start + 1
                if feature.strand == '-':
                    feature.strand = '+'
                elif feature.strand == '+':
                    feature.strand = '-'
            feature.start = new_start
            feature.end = new_end
            # if this feature is the target that the fasta file was subset around, add the word TARGET in the source column for plotting
            feature_id = feature.attributes.get('ID', [''])[0]
            if target_gene_id in feature_id:
                feature.source = 'TARGET'
            elif any([neighbor_id in feature_id for neighbor_id in neighbor_gene_ids]):
                feature.source = 'NEIGHBOR'
            else:
                feature.source = 'NONE'
            gff_out.write(str(feature) + '\n')


# def subset_fasta_and_gff(original_fasta_file, original_gff_file, scaffold, start, end, new_fasta_file, new_gff_file, target_gene_name, neighbor_gene_names, margin = 5000):
#     # subset a fasta to contain only the specified scaffold, start, and stop
#     # then, subset the gff features to contain only the features in that region
#     # update gff feature coordinates to match the new fasta
#     # first, remove the cds- prefix or the -T1.cds suffix from the gene name
#     target_gene_id = trim_cds_name(target_gene_name)
#     neighbor_gene_ids = [trim_cds_name(name) for name in neighbor_gene_names]
#     with open(new_fasta_file, 'w') as fasta_out:
#         for record in SeqIO.parse(original_fasta_file, 'fasta'):
#             if record.id == scaffold:
#                 # record the original length of the scaffold for gff coordinate updating
#                 original_length = len(record.seq)
#                 # adjust the start and end to include the margin (if possible)
#                 start = max(1,start - margin)
#                 end = min(original_length,end + margin)
#                 # subset the sequence to the start and end
#                 subset_seq = record.seq[start:end]
#                 # write the subset sequence to the new fasta file
#                 SeqIO.write(SeqIO.SeqRecord(subset_seq, id=record.id, description=''), fasta_out, 'fasta')
#     print(f'Subset {original_fasta_file} from {scaffold}:{start}-{end}')
#     with open(new_gff_file, 'w') as gff_out:
#         gff_db = gff.create_db(original_gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
#         seq_offset = start - 1
#         for feature in gff_db.region(seqid=scaffold, start=start, end=end):
#             if feature.end < start or feature.start > end:
#                 continue
#             new_start = feature.start - seq_offset
#             new_end = feature.end - seq_offset
#             new_start = max(1,new_start)
#             new_end = min(end - seq_offset,new_end)
#             if new_start >= new_end or new_end <= 1:
#                 continue
#             feature.start = new_start
#             feature.end = new_end
#             # if this feature is the target that the fasta file was subset around, add the word TARGET in the source column for plotting
#             feature_id = feature.attributes.get('ID', [''])[0]
#             if target_gene_id in feature_id:
#                 feature.source = 'TARGET'
#             elif any([neighbor_id in feature_id for neighbor_id in neighbor_gene_ids]):
#                 feature.source = 'NEIGHBOR'
#             else:
#                 feature.source = 'NONE'
#             gff_out.write(str(feature) + '\n')

def find_gf_from_gene(gene_name, isolate, pan_df):
    # find the gene family that corresponds to a given gene name in the pan_df
    for index, row in pan_df.iterrows():
        if gene_name in row[isolate]:
            return row['Gene']
    return None

def find_constant_gene_families(pan_df, b8441_gff_path):
    # for each scaffold in b8441, find a gene family that is present in all isolates and on that scaffold
    constant_gf = {}
    b8441_db = gff.create_db(b8441_gff_path,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    b8441_col = pan_df['b8441']
    for feature in b8441_db.features_of_type('CDS'):
        scaffold = feature.seqid
        if scaffold in constant_gf:
            continue
        gene_name = feature.attributes.get('ID', [''])[0]
        gene_family_row = b8441_col[b8441_col == gene_name].index[0]
        # use this value to get the corresponding row index in pan_df
        gene_family_series = pan_df.loc[gene_family_row].squeeze()
        if all([x != '' and ',' not in x and 'refound' not in x and 'pseudo' not in x for x in gene_family_series.tolist()]):
            # this should ensure the chosen gene family has no paralogs, no refound genes, and no pseudogenes
            constant_gf[scaffold] = gene_family_series
    return(constant_gf)

def master_subset(gene_families, pangenome_path, neighbor_db_path, fasta_dir, gff_dir, filter_file, margin, output_dir):
    pan = pd.read_csv(pangenome_path, sep = ',', keep_default_na=False)
    pan = pan.loc[:, pan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    ffdf = pd.read_csv(filter_file, sep = '\t', header = None)
    # convert the one-column filter file to a list
    fflist = ffdf.iloc[:,0].tolist()
    for gf in gene_families:
        isolates_with_gf, isolates_without_gf = select_isolates(gf,pan,fflist)
        isolates_to_use = isolates_with_gf + isolates_without_gf
        print(f'Selected isolates with {gf}:{",".join(isolates_with_gf)}')
        print(f'Selected isolates without {gf}:{",".join(isolates_without_gf)}')
        gf_output_dir = os.path.join(output_dir, gf)
        os.makedirs(gf_output_dir, exist_ok = True)
        constant_gf = None
        for isolate in isolates_to_use:
            neighbor_file  = os.path.join(neighbor_db_path, f'{isolate}_accessory_neighborfile.tsv')
            fasta_file = os.path.join(fasta_dir, f'{isolate}.scaffolds.fa')
            gff_file = os.path.join(gff_dir, f'{isolate}.gff3')
            # distinguish between isolates with and without the gene family
            isolate_scaffold_prefix = f'ABSENT_{isolate}'
            if isolate in isolates_with_gf:
                isolate_scaffold_prefix = f'PRESENT_{isolate}'
            if not os.path.exists(neighbor_file) or not os.path.exists(fasta_file) or not os.path.exists(gff_file):
                print(f'Could not find files for {isolate} at {neighbor_file} or {fasta_file} or {gff_file}')
                continue
            neighbor_df = pd.read_csv(neighbor_file, sep = '\t')
            target_data, neighbor1_data, neighbor2_data = find_position(gf, neighbor_df)
            positions_to_use = target_data
            if target_data[1] == 'absent':
                positions_to_use = neighbor1_data
            isolate_scaffold = positions_to_use[1]
            if isolate_scaffold == 'absent':
                print(f'Could not find scaffold for {isolate} and {gf}')
                isolates_to_use.remove(isolate)
                continue
            else:
                print(f'Found {gf} in {isolate} at {isolate_scaffold}')
            isolate_start = int(positions_to_use[2])
            isolate_end = int(positions_to_use[3])
            # subset the fasta file to include only this scaffold of interest
            scaffold_fasta_output = os.path.join(gf_output_dir, f'{isolate}_{isolate_scaffold}.fa')
            scaffold_gff_output = os.path.join(gf_output_dir, f'{isolate}_{isolate_scaffold}.gff')
            target_gene_name = target_data[0]
            neighbor_gene_names = [neighbor1_data[0], neighbor2_data[0]]
            if constant_gf is None:
                constant_gf = find_gf_from_gene(neighbor_gene_names[0], isolate, pan)
                constant_gene = neighbor_gene_names[0]
            else:
                constant_gene_str = pan[pan['Gene'] == constant_gf][isolate].values[0]
                if neighbor_gene_names[0] in constant_gene_str:
                    constant_gene = neighbor_gene_names[0]
                elif neighbor_gene_names[1] in constant_gene_str:
                    constant_gene = neighbor_gene_names[1]
                else:
                    constant_gene = None
            subset_fasta_and_gff_v2(fasta_file, gff_file, isolate_scaffold, isolate_start, isolate_end, scaffold_fasta_output, scaffold_gff_output, target_gene_name, neighbor_gene_names, constant_gene, margin)
            # in both files, replace all instances of the scaffold name with isolate_scaffold_prefix to avoid issues with shared names
            new_scaffold_name = f'{isolate_scaffold_prefix}_{isolate_scaffold}'
            subprocess.run(f'sed -i "s/{isolate_scaffold}/{new_scaffold_name}/g" {scaffold_fasta_output}', shell = True)
            subprocess.run(f'sed -i "s/{isolate_scaffold}/{new_scaffold_name}/g" {scaffold_gff_output}', shell = True) 
        # align all of the scaffolds using minimap2 and output a paf file
        scaffold_fasta_files = [os.path.join(gf_output_dir, f) for f in os.listdir(gf_output_dir) if f.endswith('.fa')]
        # perform all pairwise comparisons and concatenate the results into a single paf file
        paf_output = os.path.join(gf_output_dir, f'{gf}_scaffold_alignments.coords')
        with open(paf_output, 'w') as out:
            for i in range(len(scaffold_fasta_files)):
                for j in range(i+1, len(scaffold_fasta_files)):
                    file1 = scaffold_fasta_files[i]
                    file2 = scaffold_fasta_files[j]
                    #result = subprocess.run(['minimap2', '-x', 'asm10', '-c', file1, file2], capture_output=True, text=True)
                    # use nucmer instead of minimap2 for better handling of large indels
                    temp_delta_file = os.path.join(gf_output_dir, 'temp')
                    subprocess.run(['nucmer', '-p', temp_delta_file, file1, file2])
                    result = subprocess.run(['show-coords', '-THrcld', f'{temp_delta_file}.delta'], capture_output=True, text=True)
                    out.write(result.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--gene_families','-gf',type=str,nargs='+',
        help='''Provide a list of gene families to make files for.''',
        default=[]
        )
    parser.add_argument(
        '--gene_families_file','-gffile',type=str,
        help='''Provide a one-column csv with no header containing the names of gene families to make files for.''',
        default=None
        )
    parser.add_argument(
        '--pangenome','-p',type=str,
        help='''Provide the path to the pangenome file. This file should contain values for gene presence and absence rather than names.''',
        default=None
        )
    parser.add_argument(
        '--neighbor_db','-n',type=str,
        help='''Provide the path to the neighbor database file. This should match the clade used in the filter file.''',
        default=None
        )
    parser.add_argument(
        '--fasta_dir','-fd',type=str,
        help='''Provide the path to the fasta directory.''',
        default=None
        )
    parser.add_argument(
        '--gff_dir','-gd',type=str,
        help='''Provide the path to the gff directory.''',
        default=None
        )
    parser.add_argument(
        '--filter_file','-ff',type=str,
        help='''Provide the path to the filter file. This should be a one-column file of isolate names from the same clade used in the neighbor database.''',
        default=None
        )
    parser.add_argument(
        '--margin','-m',type=int,
        help='''Provide the margin around the target gene.''',
        default=8000
        )
    parser.add_argument(
        '--output_dir','-o',type=str,
        help='''Provide the path to the output directory.''',
        default=None
        )
    args = parser.parse_args()
    for path in [args.pangenome, args.neighbor_db, args.fasta_dir, args.gff_dir, args.filter_file]:
        if not os.path.exists(path):
            print(f'Could not find file at {path}')
            quit(1)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if args.gene_families_file is not None:
        if not os.path.exists(args.gene_families_file):
            print(f'Could not find gene families file at {args.gene_families_file}')
            quit(1)
        else:
            gf_df = pd.read_csv(args.gene_families_file, header = None)
            args.gene_families.extend(gf_df.iloc[:,0].tolist())
    master_subset(args.gene_families, args.pangenome, args.neighbor_db, args.fasta_dir, args.gff_dir, args.filter_file, args.margin, args.output_dir)

if __name__ == '__main__':
    main()
