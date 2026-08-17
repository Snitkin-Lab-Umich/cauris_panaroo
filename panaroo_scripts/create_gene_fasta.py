import os
import subprocess
from Bio import SeqIO
import argparse
import gffutils as gff
import pandas as pd

def extract_sequences(input_file, assemblies_dir, gff_dir, pangenome_file, output_file):
    # take an input file as a two-column csv, with 'gene_family' containing the family and 'clades' containing a string of the clade data
    # use the pangenome file to find the gene names for each family, based on the provided isolate name
    # use the gff file to find the location of each gene
    # use the assembly fasta file to extract the sequence of each gene
    # write all of these genes to a single output fasta file
    genefam_data = pd.read_csv(input_file)
    fullpan = pd.read_csv(pangenome_file,keep_default_na=False)
    pan = fullpan.loc[:, fullpan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    with open(output_file, 'w') as out_fasta:
        for index, row in genefam_data.iterrows():
            genefam = row['gene_family']
            clade_data = row['clades']
            isolate_gene_dict = make_isolate_dict(genefam, pan)
            if isolate_gene_dict == {}:
                print(f'Warning: no searchable genes found in gene family {genefam} at all')
                continue
            isolate, gene_name = find_best_gene(isolate_gene_dict, gff_dir)
            # remove the .cds suffix from the gene name to extract the mRNA feature rather than the CDS feature
            # if using the spliced data, use the CDS name instead
            #gene_name_mrna = gene_name.replace('.cds', '')
            gene_name_mrna = gene_name
            gff_file = os.path.join(gff_dir, f'{isolate}.gff3')
            if not os.path.isfile(gff_file):
                gff_file = os.path.join(gff_dir, f'{isolate}.gff')
            if not os.path.isfile(gff_file):
                print(f'Unable to locate gff file for isolate {isolate} at {gff_file}')
                continue
            fasta_file = os.path.join(assemblies_dir, f'{isolate}.scaffolds.fa')
            query_fasta_record = extract_record_from_gff(gff_file, fasta_file, gene_name_mrna)
            #print(f'Found gene {query_fasta_record.id} in isolate {isolate} for gene family {genefam}')
            #print(f'gff: {gff_file}, fasta: {fasta_file}, gene_name: {gene_name_mrna}')
            # write the record to the output fasta file
            query_fasta_record.id = f'{isolate}_{genefam}_{gene_name_mrna}_{clade_data}'
            query_fasta_record.description = ''
            SeqIO.write(query_fasta_record, out_fasta, 'fasta')

def make_isolate_dict(genefam, pan):
    if genefam not in pan['Gene'].values:
        print(f'{genefam} not found in pangenome file')
        return {}
    genefam_row_index = pan[pan['Gene'] == genefam].index[0]
    pan_data_series = pan.loc[genefam_row_index]
    # create a dictionary of isolate:gene_name for the isolates that have this gene family
    isolate_gene_dict = {}
    for isolate in pan_data_series.index:
        if isolate == 'Gene':
            continue
        gene_name = pan_data_series[isolate]
        if gene_name == '' or 'refound' in gene_name or 'pseudo' in gene_name:
            continue
        # method to keep pseudo genes in the search
        # if gene_name.endswith('_pseudo'):
        #     gene_name = gene_name.replace('_pseudo','')
        isolate_gene_dict[isolate] = gene_name
    return isolate_gene_dict

def extract_record_from_gff(gff_file, fasta_file, gene_name):
    # parse the gff file to get the start and end positions of the gene
    db = gff.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    gene_feature = db[gene_name]
    seqid = gene_feature.seqid
    start = gene_feature.start
    end = gene_feature.end
    strand = gene_feature.strand
    # extract the sequence from the fasta file
    #print(f'Extracting sequence for gene {gene_name} from {fasta_file}, seqid: {seqid}, start: {start}, end: {end}, strand: {strand}')
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == seqid:
            gene_seq = record.seq[start-1:end]  # gff is 1-based, python is 0-based
            gene_record = SeqIO.SeqRecord(gene_seq, id=gene_name, description='')
            return(gene_record)
    print(f'Could not find sequence for gene {gene_name} in fasta file {fasta_file}')
    quit(1)



def find_best_gene(isolate_gene_dict, gff_dir):
    # take a dict of isolate:gene_name
    # eventually, this will need to compare lengths/quality of sequences
    # for now, simply prioritize any isolates that start with UM_ or Chi_, as these are hybrid assemblies
    final_isolate = None
    final_gene_name = None
    for isolate in isolate_gene_dict:
        if isolate.startswith('UM_') or isolate.startswith('Chi_'):
            final_isolate = isolate
            final_gene_name = isolate_gene_dict[isolate]
            break
        final_isolate = isolate
        final_gene_name = isolate_gene_dict[isolate]
    # some gene names are actually lists of multiple genes separated by a semicolon
    if ';' in final_gene_name:
        final_gene_name = final_gene_name.split(';')[0]
    return final_isolate, final_gene_name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a one-column csv containing the names of accessory gene families to extract sequences for. These should be the exact group names from Panaroo.''',
        default=None
        )
    parser.add_argument(
        '--assemblies','-a',type=str,
        help='''Provide the name of a directory containing assemblies for all isolates. This should be in the same format as the Panaroo input, with
        file names ending in .scaffolds.fa''',
        default=None
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide the name of a directory containing gff files for all isolates. This should be in the same format as the Panaroo input, with
        file names ending in .gff3''',
        default=None
        )
    parser.add_argument(
        '--pangenome','-p',type=str,
        help='''Provide a path to the presence/absence matrix provided by Panaroo. Gene names should be present in this file.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output file. This will be a fasta file with a single entry for each gene family in the input list.''',
        default=None
        )
    args = parser.parse_args()
    extract_sequences(args.input, args.assemblies, args.gff, args.pangenome, args.output)


if __name__ == '__main__':
    main()