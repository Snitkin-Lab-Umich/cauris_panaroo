import os
import argparse
import gffutils as gff
import pandas as pd

def get_all_gff_annotation_terms(gff_dir, pangenome_pa, annotation_types, output_file, skip_isolates=['b8441']):
    # start with a directory of gff files and a pangenome presence absence file with gene families as rows and genomes as columns
    # for each isolate, extract all annotations of the specified types, then add these annotations to the corresponding gene families
    # repeat until no gene families are left without annotations
    # output a table with one row for each gene family, with each of the specified annotation types as an additional column
    fullpan = pd.read_csv(pangenome_pa, sep = ',', keep_default_na=False)
    pan = fullpan.loc[:, fullpan.columns.str.startswith(("Gene", "SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441"))]
    annotation_dict = {'GeneFamily':[], **{annotation_type:[] for annotation_type in annotation_types}, 'AnnotationSource':[], 'GeneName':[], 'Scaffold':[],'Start':[],'End':[]}
    print(f'columns in annotation dict: {[','.join(annotation_dict.keys())]}')
    genes_to_annotate = pan['Gene'].tolist()
    isolate_list = pan.columns[1:]
    # order this list with b8441 first, then UM_Caur_4, then Chi_Caur_3, then any remaining UM or Chi isolates, then the rest in any order
    isolate_list = sorted(isolate_list, key=lambda x: (0 if x == 'b8441' else 1, 0 if x == 'UM_Caur_4' else 1, 0 if x == 'Chi_Caur_3' else 1, 0 if x.startswith('UM') or x.startswith('Chi') else 1, x))
    for skip_isolate in skip_isolates:
        if skip_isolate in isolate_list:
            isolate_list.remove(skip_isolate)
    isolate_index = 0
    while isolate_index < len(isolate_list) and len(genes_to_annotate) > 0:
        isolate = isolate_list[isolate_index]
        isolate_index += 1
        print(f'Starting isolate {isolate}')
        print(f'{len(genes_to_annotate)} gene families left to annotate.')
        print(f'length of each column in annotation dict: {','.join([str(len(annotation_dict[x])) for x in annotation_dict.keys()])}')
        gff_file = os.path.join(gff_dir, f'{isolate}.gff3')
        if not os.path.exists(gff_file):
            print(f'Could not find gff file for {isolate} at {gff_file}')
            continue
        gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
        for feature in gff_db.features_of_type('CDS'):
            #print(f'Processing {feature.id} from {isolate}')
            cds_name = feature.attributes.get('ID', ['NA'])[0]
            # get the parent of this CDS, which should be the mRNA name, and use it to extract the mRNA feature as well
            parent_name = feature.attributes.get('Parent', ['NA'])[0]
            mrna_feature = gff_db[parent_name]
            annotation_list = mrna_feature.attributes.get('Dbxref', ['']) + mrna_feature.attributes.get('Note', ['']) + mrna_feature.attributes.get('note', [''])
            # determine which gene family contains this gene_name, if any
            gene_family = pan['Gene'].loc[pan[isolate].str.contains(cds_name)].tolist()
            # if this gene family has exactly one match and is still present in genes_to_annotate, add it to the annotation dict
            if len(gene_family) == 1 and gene_family[0] in genes_to_annotate:
                gene_family = gene_family[0]
                #print(f'Found gene family {gene_family} for {cds_name} from {isolate}')
                annotation_dict['GeneFamily'].append(gene_family)
                annotation_dict['AnnotationSource'].append(isolate)
                annotation_dict['GeneName'].append(cds_name)
                annotation_dict['Scaffold'].append(feature.seqid)
                annotation_dict['Start'].append(feature.start)
                annotation_dict['End'].append(feature.end)
                for annotation_type in annotation_types:
                    to_append = []
                    for ann_entry in annotation_list:
                        if annotation_type in ann_entry:
                            to_append.append(ann_entry)
                    if to_append == []:
                        annotation_dict[annotation_type].append('NA')
                    else:
                        annotation_dict[annotation_type].append(','.join(to_append))
                genes_to_annotate.remove(gene_family)
    annotation_df = pd.DataFrame(annotation_dict)
    annotation_df.to_csv(output_file, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--gff_dir','-g',type=str,
        help='''Provide the path to the input gff directory.''',
        default=None
        )
    parser.add_argument(
        '--pangenome','-p',type=str,
        help='''Provide the path to the pangenome file. This file should contain names of genes rather than numbers.''',
        default=None
        )
    parser.add_argument(
        '--annotation_types','-at',type=str,nargs='+',
        help='''Provide a list of annotation types to include. Make sure this matches the format of the annotations in the gff file.''',
        default=['InterPro', 'PFAM', 'COG', 'EggNog', 'BUSCO']
        )
    parser.add_argument(
        '--skip_isolates','-si',type=str,nargs='+',
        help='''Provide a list of isolates to skip.''',
        default=['b8441']
        )
    parser.add_argument(
        '--output_file','-o',type=str,
        help='''Provide the path to the output file.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.exists(args.gff_dir):
        print(f'Could not find gff directory at {args.gff_dir}')
        quit(1)
    if not os.path.exists(args.pangenome):
        print(f'Could not find pangenome file at {args.pangenome}')
        quit(1)
    get_all_gff_annotation_terms(args.gff_dir, args.pangenome, args.annotation_types, args.output_file, args.skip_isolates)

if __name__ == '__main__':
    main()