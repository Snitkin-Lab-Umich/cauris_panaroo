import os
import argparse
import gffutils as gff
from Bio import SeqIO
import pandas as pd

def add_annotations(input_summary_file, gff_dir_original, pangenome_file, assemblies_dir, output_file):
    # steps for each gene family
    # 1) read the pangenome matrix, get the row that matches the gene family name
    # 2) get the list of isolates that actually have that gene family
    # 3) for each isolate, parse the gff file to get the start and end positions of that gene
    # 4) extract this sequence from the corresponding fasta file
    # 5) compare the extracted sequences and write the longest/best one to the output file
    #####
    genefam_data = pd.read_csv(input_summary_file, sep='\t')
    fullpan = pd.read_csv(pangenome_file, sep = ',', keep_default_na=False)
    additional_genefam_data = {'annotation_source': [],'gene_id': [], 'name': [], 'product': [], 'annotation_data': [], 'sequence': []}
    # subset pan to only contain columns for the gene family names and isolates
    fullpan = fullpan.loc[:, fullpan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    for genefam in genefam_data['GeneFamily']:
        isolate_gene_dict = make_isolate_dict(genefam, fullpan)
        if isolate_gene_dict == {}:
            print(f'Error: no searchable genes found in gene family {genefam}')
            quit(1)
        isolate, gene_name = find_best_gene(isolate_gene_dict)
        # remove the .cds suffix from the gene name to extract the mRNA feature rather than the CDS feature
        # if using the spliced data, use the CDS name instead
        #gene_name_mrna = gene_name.replace('.cds', '')
        gene_name_mrna = gene_name
        gff_file_original = os.path.join(gff_dir_original, f'{isolate}.gff3')
        if not os.path.isfile(gff_file_original):
            gff_file_original = os.path.join(gff_dir_original, f'{isolate}.gff')
        if not os.path.isfile(gff_file_original):
            print(f'Unable to locate gff files for isolate {isolate}')
            continue
        # if an assebly directory was provided, extract the sequence as well as the annotation data
        if assemblies_dir is not None:
            fasta_file = os.path.join(assemblies_dir, f'{isolate}.scaffolds.fa')
            gene_id, name, product, annotation_data, gene_record = extract_annotation_from_gff(gff_file_original, gene_name_mrna, fasta_file)
            additional_genefam_data['sequence'].append(str(gene_record.seq))
        else:
            gene_id, name, product, annotation_data = extract_annotation_from_gff(gff_file_original, gene_name_mrna)
        # add all annotation data to the additional_genefam_data dictionary
        additional_genefam_data['annotation_source'].append(isolate)
        additional_genefam_data['gene_id'].append(gene_id)
        additional_genefam_data['name'].append(name)
        additional_genefam_data['product'].append(product)
        additional_genefam_data['annotation_data'].append(annotation_data)
    # add the additional annotation data to the original summary file and write to a new output file
    genefam_data['annotation_source'] = additional_genefam_data['annotation_source']
    genefam_data['gene_id'] = additional_genefam_data['gene_id']
    genefam_data['name'] = additional_genefam_data['name']
    genefam_data['product'] = additional_genefam_data['product']
    genefam_data['annotation_data'] = additional_genefam_data['annotation_data']
    if assemblies_dir is not None:
        genefam_data['sequence'] = additional_genefam_data['sequence']
    genefam_data.to_csv(output_file, index=False, sep='\t')


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


def find_best_gene(isolate_gene_dict):
    # take a dict of isolate:gene_name
    # eventually, this will need to compare lengths/quality of sequences
    # for now, simply prioritize any isolates that start with UM_ or Chi_, as these are hybrid assemblies
    # also, give even higher priority to the three reference isolates (b8441, UM_Caur_4, Chi_Caur_3)
    final_isolate = None
    final_gene_name = None
    for isolate in isolate_gene_dict:
        if isolate in ['b8441','UM_Caur_4','Chi_Caur_3']:
            final_isolate = isolate
            final_gene_name = isolate_gene_dict[isolate]
            break
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


def extract_annotation_from_gff(gff_file, gene_name, fasta_file=None):
    # parse the gff file to get the annotation data for the gene
    # note that the gene name is the name of the CDS feature, but the annotation data is actually stored in the parent mRNA and parent gene features
    db = gff.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    gene_feature_cds = db[gene_name]
    gene_feature_mrna = db[gene_feature_cds.attributes['Parent'][0]]
    gene_feature_gene = db[gene_feature_mrna.attributes['Parent'][0]]
    product = gene_feature_mrna.attributes.get('product', ['NA'])[0]
    dbxref = ','.join(gene_feature_mrna.attributes.get('Dbxref', ['NA']))
    note = ','.join(gene_feature_mrna.attributes.get('note', ['NA']))
    gene_id = gene_feature_gene.attributes.get('ID', ['NA'])[0]
    name = gene_feature_gene.attributes.get('Name', ['NA'])[0]
    # combine dbxref and note into a single string with all annotation data
    annotation_data = f'{dbxref},{note}'
    if fasta_file is None:
        return gene_id, name, product, annotation_data
    else:
        start = gene_feature_mrna.start
        end = gene_feature_mrna.end
        seqid = gene_feature_mrna.seqid
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id == seqid:
                gene_seq = record.seq[start-1:end]  # gff is 1-based, python is 0-based
                gene_record = SeqIO.SeqRecord(gene_seq, id=gene_name, description='')
                return gene_id, name, product, annotation_data, gene_record



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a summary file of validated accessory genes. This file should contain a column named "GeneFamily" that matches the gene family names in the pangenome file.''',
        default=None
        )
    parser.add_argument(
        '--gff_dir','-g',type=str,
        help='''Provide the path to the directory containing the GFF files.''',
        default=None
        )
    parser.add_argument(
        '--assemblies_dir','-a',type=str,
        help='''(Optional) Provide the path to the directory containing the assembly fasta files.''',
        default=None
        )
    parser.add_argument(
        '--pangenome','-p',type=str,
        help='''Provide a path to the presence/absence matrix provided by Panaroo. Gene names should be present in this file.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''(Optional) Provide the path to the output file.''',
        default='validated_accessory_gene_summary.tsv'
        )
    args = parser.parse_args()
    if not args.gff_dir.endswith('/'):
        args.gff_dir += '/'
    if not os.path.exists(args.gff_dir):
        print(f'Cannot locate GFF directory at {args.gff_dir}.')
        quit(1)
    if args.assemblies_dir is not None:
        if not args.assemblies_dir.endswith('/'):
            args.assemblies_dir += '/'
        if not os.path.exists(args.assemblies_dir):
            print(f'Cannot locate assemblies directory at {args.assemblies_dir}.')
            quit(1)
    add_annotations(args.input, args.gff_dir, args.pangenome, args.assemblies_dir, args.output)

if __name__ == '__main__':
    main()