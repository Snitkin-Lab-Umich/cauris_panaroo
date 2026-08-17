import os
from Bio import SeqIO
import argparse
import gffutils as gff
import pandas as pd

def read_pangenome_file(pangenome_file, gene_family):
    # take a pangenome dataframe with the names of genes
    # return a dict of isolate:gene_name for the given gene family
    pan = pd.read_csv(pangenome_file,keep_default_na=False)
    pan = pan.loc[:, pan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    outdict = {}
    if gene_family not in pan['Gene'].values:
        print(f'{gene_family} not found in pangenome file')
        quit(1)
    gene_family_row_index = pan[pan['Gene'] == gene_family].index[0]
    pan_data_series = pan.loc[gene_family_row_index]
    for isolate in pan_data_series.index:
        if isolate == 'Gene':
            continue
        gene_name = pan_data_series[isolate]
        if gene_name == '' or 'refound' in gene_name or 'pseudo' in gene_name:
            continue
        outdict[isolate] = gene_name
    return outdict

def read_master_qc_file(master_qc_file):
    m = pd.read_csv(master_qc_file, keep_default_na=False)
    outdict = {}
    for index, row in m.iterrows():
        isolate_name = row['Sample']
        clade_name = row['auriclass_clade']
        clade_name = clade_name.replace('Clade ', 'Clade_')
        outdict[isolate_name] = clade_name
    return outdict

def make_protein_fasta(gene_family, pangenome_file, fasta_dir, master_qc, output_file):
    # take a gene family name, a pangenome file, and a directory of protein fasta files
    # use these to create a multi-fasta file with one entry for each isolate that has a gene in this family
    isolate_dict = read_pangenome_file(pangenome_file, gene_family)
    clade_dict = read_master_qc_file(master_qc)
    with open(output_file, 'w') as out_fasta:
        for isolate_name in isolate_dict:
            gene_name = isolate_dict[isolate_name]
            protein_name = gene_name.replace('.cds', '')
            clade_name = clade_dict.get(isolate_name, 'Unknown_Clade')
            protein_fasta = os.path.join(fasta_dir,f'{isolate_name}.proteins.fa')
            if not os.path.isfile(protein_fasta):
                print(f'Could not find protein fasta file for isolate {isolate_name} at {protein_fasta}')
                continue
            for record in SeqIO.parse(protein_fasta, 'fasta'):
                if record.id == protein_name:
                    record.id = f'{isolate_name}_{clade_name}_{protein_name}'
                    record.description = ''
                    SeqIO.write(record, out_fasta, 'fasta')
                    break

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--pangenome_file','-p',type=str,
        help='''Provide the path to a pangenome file (with gene names).''',
        default=None
        )
    parser.add_argument(
        '--fasta_dir','-fd',type=str,
        help='''Provide the path to a directory of protein fasta files, with file names ending in .proteins.fa''',
        default=None
        )
    parser.add_argument(
        '--gene_family','-gf',type=str,
        help='''Provide the name of the gene family to extract sequences for.''',
        default=None
        )
    parser.add_argument(
        '--master_qc','-mq',type=str,
        help='''Provide a path to a master QC file with clade names.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output file.''',
        default=None
        )
    args = parser.parse_args()
    make_protein_fasta(args.gene_family, args.pangenome_file, args.fasta_dir, args.master_qc, args.output)


if __name__ == '__main__':
    main()