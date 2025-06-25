import os
import argparse
import gffutils as gff

# parse a gff3 file, looking for matches for a list of genes
# this will ONLY work for the gff files from funannotate - the b8441 reference has a different format
def get_gff_genes(gff_file,genelist,isolatename):
    outdata = []
    gffdb = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    for gene in genelist:
        genedata = gff_finder(gffdb,gene)
        if genedata == 'NA':
            print(f'Gene {gene} not found in {gff_file}')
            continue
        if genedata.featuretype == 'gene':
            print(f'Gene feature provided at {gene} instead of an mRNA feature.')
        elif genedata.featuretype != 'mRNA':
            while genedata.featuretype != 'mRNA':
                genedata = gffdb[genedata.attributes['Parent'][0]]
                # this will go up the hierarchy until it finds the parent mRNA feature
        product = gff_finder(genedata.attributes, 'product')
        if product != 'NA':
            product = product[0]
        dbxref = gff_finder(genedata.attributes, 'Dbxref')
        if dbxref != 'NA':
            dbxref = ','.join(dbxref)
        outdata.append([isolatename, gene, genedata.attributes['ID'][0], min(genedata.start,genedata.end), max(genedata.start,genedata.end), product, dbxref])
    return outdata

def get_gff_genes_b8441(gff_file,genelist,isolatename):
    outdata = []
    gffdb = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    for gene in genelist:
        genedata = gff_finder(gffdb,gene)
        if genedata == 'NA':
            print(f'Gene {gene} not found in {gff_file}')
            continue
        if genedata.featuretype != 'CDS':
            print(f'CDS feature not provided at {gene} for B8441!')
        product = gff_finder(genedata.attributes, 'product')
        if product != 'NA':
            product = product[0]
        dbxref = gff_finder(genedata.attributes, 'Dbxref')
        if dbxref != 'NA':
            dbxref = ','.join(dbxref)
        outdata.append([isolatename, gene, genedata.attributes['ID'][0], min(genedata.start,genedata.end), max(genedata.start,genedata.end), product, dbxref])
    return outdata

def gff_finder(gff_feature, attribute_name):
    # extracts the value of a specific attribute from a GFF feature
    try:
        return gff_feature[attribute_name]
    except gff.FeatureNotFoundError:
        return 'NA'
    except KeyError:
        return 'NA'


def read_all_gff(input_file,output_file,gff_dir):
    with open(input_file, 'r') as fhin, open(output_file, 'w') as fhout:
        gff_dict = {}
        next(fhin)
        for line in fhin:
            gff_file, gene = line.strip().split('\t')
            if ';' in gene:
                gene = gene.split(';')
            else:
                gene = [gene]
            if gff_file not in gff_dict:
                gff_dict[gff_file] = gene
            else:
                gff_dict[gff_file] += gene
        fhout.write('\t'.join(['Isolate', 'Gene', 'mRNA_ID', 'Start', 'End', 'product', 'dbxref']) + '\n')
        for gff_file in gff_dict:
            isolatename = gff_file.split('/')[-1].split('.gff')[0]
            gff_file_path = gff_dir + gff_file
            if not os.path.exists(gff_file_path):
                print(f'GFF file not found at {gff_file_path}.')
                quit(1)
            if 'b8441' in isolatename or 'B8441' in isolatename:
                writedata = get_gff_genes_b8441(gff_file_path, gff_dict[gff_file], isolatename)
            else:
                writedata = get_gff_genes(gff_file_path, gff_dict[gff_file], isolatename)
            for writeline in writedata:
                fhout.write('\t'.join([str(x) for x in writeline]) + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a tab-separated file, with the first column being the isolate name and the second column being the feature name.
        mRNA feature names are preferred.''',
        default=None
        )
    parser.add_argument(
        '--gff_dir','-g',type=str,
        help='''Provide the path to the directory containing the GFF files.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''(Optional) Provide the path to the output file.''',
        default='gene_locations.tsv'
        )
    args = parser.parse_args()
    if not args.gff_dir.endswith('/'):
        args.gff_dir += '/'
    if not os.path.exists(args.gff_dir):
        print(f'Cannot locate GFF directory at {args.gff_dir}.')
        quit(1)
    read_all_gff(args.input,args.output,args.gff_dir)


if __name__ == '__main__':
    main()