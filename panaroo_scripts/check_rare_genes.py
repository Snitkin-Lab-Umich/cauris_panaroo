import os
import argparse
import networkx as nx
import gffutils as gff
import pandas as pd

def find_adjacent_genes(gene_families, gml_file):
    # take a pangenome graph in gml format and list of gene families
    # return a dictionary of gene families and their adjacent genes
    pangml = nx.read_gml(gml_file, label = 'name')
    adjacent_genes = {}
    for family in gene_families:
        if family in pangml:
            neighbors = list(pangml.neighbors(family))
            adjacent_genes[family] = neighbors
        else:
            adjacent_genes[family] = []
    return adjacent_genes

# make sure to check that this returns the expected structure and cds format
def parse_gene_families(family_file, isolate_prefixes):
    # take a file where the first column is the gene family and the header is the isolate
    # return a dictionary of gene_family:isolate:family_names
    # this is intended for the full table of all gene families, not the table of rare genes
    gene_families = {}
    with open(family_file, 'r') as fh:
        # the header contains the isolate names, as well as which columns to skip
        header = fh.readline().strip().split(',')
        for line in fh:
            lineparts = line.strip().split(',')
            family_name = lineparts[0]
            gene_families[family_name] = {}
            for i in range(1, len(lineparts)):
                isolate = header[i]
                if isolate.startswith(isolate_prefixes):
                    gene_names = lineparts[i]
                    gene_families[family_name][isolate] = []
                    # note that gene_names may be multiple cds sequences, a refound gene, a pseudogene, or absent
                    # first, split the gene names in the case that there are multiple
                    gene_names_list = gene_names.split(';')
                    for gene in gene_names_list:
                        # check if this is a refound gene
                        if 'refound' in gene and 'cds' not in gene:
                            gene_families[family_name][isolate].append('refound')
                        elif 'pseudo' in gene:
                            # it is possible to have refound pseudogenes or existing gene marked as pseudogenes
                            # refound pseudogenes should be removed by the lines above, so everything here should have the format of:
                            # FUN_005627-T1.cds_pseudo
                            gene_families[family_name][isolate].append(gene.replace('_pseudo', ''))
                        elif 'cds' in gene:
                            gene_families[family_name][isolate].append(gene)
                        elif gene == 'NA' or gene == '':
                            gene_families[family_name][isolate].append('absent')
                        else:
                            print(f'Unrecognized gene format: {gene}')
                            quit(1)
    return(gene_families)


# remember that there are several outcomes to check for:
# 1) the rare gene is present in the isolate, but:
#   a) one or both of the adjacent genes are absent entirely
#   b) one or both of the adjacent genes are present on a different scaffold as the rare gene
#   c) both adjacent genes are present on the same scaffold, with one upstream and one downstream
#   d) there is only one adjacent gene, and I can only check for presence and shared scaffold
#   e) there are no adjacent genes, and I cannot check for anything
#   f) there are multiple adjacent genes, and I need to find one upstream and one downstream (if this is rare, I want to skip this case)
# 2) the rare gene is absent in the isolate, and:
#   a) one or both of the adjacent genes are also absent
#   b) both of the adjacent genes are present, but one of them is on a different scaffold
#       b1) note that the location of the rare gene cannot be directly checked in this specific case
#   c) both adjacent genes are present on the same scaffold, with one upstream and one downstream
#   d) there is only one adjacent gene, and I can only check for presence
#   e) there are no adjacent genes, and I cannot check for anything
#   f) there are multiple adjacent genes, and I need to find one upstream and one downstream (if this is rare, I want to skip this case)

def gff_locator(gene_family_name, gff_db, gene_families_dict, isolate_name):
    # take a gene family name, a gffutils database, a gene families dictionary, and an isolate name
    # return the scaffold, start, and end of each coding sequence in the gene family for the given isolate
    # first, convert the gene family name to a list of coding sequences
    location_dict = {'chromosome':[], 'start':[], 'end':[]}
    gene_list = gene_families_dict[gene_family_name][isolate_name]
    for gene in gene_list:
        if gene == 'absent' or gene == 'refound':
            location_dict['chromosome'].append('not_found')
            location_dict['start'].append(-1)
            location_dict['end'].append(-1)
        else:
            try:
                gene_entry = gff_db[gene]
                location_dict['chromosome'].append(gene_entry.seqid)
                location_dict['start'].append(gene_entry.start)
                location_dict['end'].append(gene_entry.end)
            except gff.exceptions.FeatureNotFoundError:
                location_dict['chromosome'].append('not_found')
                location_dict['start'].append(-1)
                location_dict['end'].append(-1)
    return location_dict

def absent_neighbor_solver(neighbor_locations):
    # take a list two neighboring gene families, assuming that the rare gene is not present
    # account for the fact that each neighboring gene family may have an arbitrary number of coding sequences
    neighbor_1 = neighbor_locations[0]
    neighbor_2 = neighbor_locations[1]
    # if either neighbor is absent or refound, no locations can be compared
    if neighbor_1['chromosome'] == ['not_found'] or neighbor_2['chromosome'] == ['not_found']:
        return 'GeneAbsent_NeighborsAbsent'
    # check if there are any shared scaffolds among the coding sequences of each neighbor
    shared_scaffolds = set(neighbor_1['chromosome']).intersection(set(neighbor_2['chromosome']))
    if len(shared_scaffolds) == 0:
        return 'GeneAbsent_NeighborsPresent_DifferentScaffold'
    else:
        # if a scaffold is shared, there is nothing more to check, since upstream and downstream cannot be determined
        return 'GeneAbsent_NeighborsPresent_SameScaffold'

def present_neighbor_solver(neighbor_locations, rare_gene_location):
    # take a list two neighboring gene families, assuming that the rare gene is present in the genome
    # account for the fact that each neighboring gene family and the rare gene itself may have an arbitrary number of coding sequences
    # note that refound genes will currently be treated the same as absent genes!
    neighbor_1 = neighbor_locations[0]
    neighbor_2 = neighbor_locations[1]
    # if there is a coding sequence in either neighbor that shares a scaffold with the rare gene, save it
    # if there are multiple, just take the first one
    # if one is found, determine if it is upstream or downstream of the rare gene
    neighbor_1_status = 'DifferentScaffold'
    if rare_gene_location['chromosome'] in neighbor_1['chromosome']:
        for i in range(len(neighbor_1['chromosome'])):
            if neighbor_1['chromosome'][i] == rare_gene_location['chromosome']:
                if neighbor_1['start'][i] < rare_gene_location['start']:
                    neighbor_1_status = 'Upstream'
                if neighbor_1['start'][i] >= rare_gene_location['start']:
                    neighbor_1_status = 'Downstream'
                break
    if neighbor_1['chromosome'] == ['not_found']:
        neighbor_1_status = 'NotPresent'
    neighbor_2_status = 'DifferentScaffold'
    if rare_gene_location['chromosome'] in neighbor_2['chromosome']:
        for i in range(len(neighbor_2['chromosome'])):
            if neighbor_2['chromosome'][i] == rare_gene_location['chromosome']:
                if neighbor_2['start'][i] < rare_gene_location['start']:
                    neighbor_2_status = 'Upstream'
                if neighbor_2['start'][i] >= rare_gene_location['start']:
                    neighbor_2_status = 'Downstream'
                break
    if neighbor_2['chromosome'] == ['not_found']:
        neighbor_2_status = 'NotPresent'
    # return the status of both neighbors
    status_string = f'{neighbor_1_status}{neighbor_2_status}'
    if status_string in ['UpstreamDownstream', 'DownstreamUpstream']:
        status_string = 'PresentFlanking'
    if status_string in ['UpstreamUpstream', 'DownstreamDownstream']:
        status_string = 'PresentNonflanking'
    return(f'GenePresent_Neighbors{status_string}')

def check_adjacent_genes(adjacent_genes, gene_families, gff_file, rare_gf_list, isolate_name):
    # adjacent_genes is a dictionary with the format gene_family:[adjacent_gene1, adjacent_gene2, etc]
    # gene_families is a dictionary with the format gene_family:isolate:[gene1, gene2, etc]
    # gff_file is a path to the gff file for a specific isolate, which contains the location of each gene
    # for each gene family in rare_gf_list, return a summary of the status of the adjacent genes
    # the output should be in the format gene_family:status
    # note that refound genes are currently treated though they are absent!
    out_dict = {}
    gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    for rgf in rare_gf_list:
        # get the list of coding sequences for this rare gene family in this isolate
        # for the sake of simplicity, only take the first coding sequence found for the rare gene
        # note that there are currently only 5 cases where a rare gene has multiple coding sequences in an isolate
        rare_gene_cds = gene_families[rgf][isolate_name][0]
        # check if the rare gene is present in this isolate
        if rare_gene_cds == 'absent' or rare_gene_cds == 'refound':
            rare_gene_status = 'absent'
        else:
            rare_gene_status = 'present'
        # get the list of adjacent genes for this rare gene family
        adj_genes = adjacent_genes.get(rgf, [])
        if len(adj_genes) < 2:
            out_dict[rgf] = f'Gene{rare_gene_status}_TooFewNeighbors'
        elif len(adj_genes) > 2:
            out_dict[rgf] = f'Gene{rare_gene_status}_TooManyNeighbors'
        else:
            # convert each neighbor to a dictionary containing the locations of all coding sequences that gene family contains
            neighbor_locations = [gff_locator(x, gff_db, gene_families, isolate_name) for x in adj_genes]
            if rare_gene_status == 'absent':
                # if the gene is absent, I can only determine if the neighboring genes are on the same scaffold
                neighbor_status = absent_neighbor_solver(neighbor_locations)
            if rare_gene_status == 'present':
                # if the gene is present, I also need its position
                try:
                    rare_gene_gff_data = gff_db[rare_gene_cds]
                    rare_gene_location = {'chromosome': rare_gene_gff_data.seqid, 'start': rare_gene_gff_data.start, 'end': rare_gene_gff_data.end}
                except gff.exceptions.FeatureNotFoundError:
                    print(f'Could not find {rare_gene_cds} in the gff file for {isolate_name}.')
                    quit(1)
                neighbor_status = present_neighbor_solver(neighbor_locations, rare_gene_location)
            out_dict[rgf] = neighbor_status
    return out_dict

def parse_rare_gene_families(rare_gene_file):
    # take a one-column file of rare gene families and convert it into a list
    rare_genes = []
    with open(rare_gene_file, 'r') as fh:
        next(fh)
        for line in fh:
            rare_genes.append(line.strip())
    return rare_genes

def get_isolate_list(family_file, isolate_prefixes):
    isolate_list = []
    with open(family_file, 'r') as fh:
        # the header contains the isolate names, as well as which columns to skip
        header = fh.readline().strip().split(',')
        for isolate in header:
            if isolate.startswith(isolate_prefixes):
                isolate_list.append(isolate)
    return isolate_list

def check_all_rare_genes(gff_dir, gene_family_file, rare_gene_file, gml_file, output_file):
    # list out the prefixes for the isolate names
    isolate_prefixes = ("SRR", "ARR", "DRR", "ERR", "UM_", "Chi_", "b8441")
    # get gene families dictionary
    gene_families = parse_gene_families(gene_family_file, isolate_prefixes)
    # get isolate list
    isolate_list = get_isolate_list(gene_family_file, isolate_prefixes)
    # get rare gene list
    rare_genes_list = parse_rare_gene_families(rare_gene_file)
    # get adjacent genes dictionary
    adjacent_genes = find_adjacent_genes(rare_genes_list, gml_file)
    # for each isolate, check each rare gene
    dfdata = {'Gene':rare_genes_list}
    for isolate in isolate_list:
        gff_file = os.path.join(gff_dir, f'{isolate}.gff')
        if not os.path.isfile(gff_file):
            print(f'Could not find gff file for {isolate} at {gff_file}.')
            quit(1)
        isolate_results = check_adjacent_genes(adjacent_genes, gene_families, gff_file, rare_genes_list, isolate)
        dfdata[isolate] = [isolate_results[rgf] for rgf in rare_genes_list]
    outdf = pd.DataFrame(dfdata)
    outdf.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--locations','-l',type=str,
        help='''Provide a path to a directory of Prokka gff files.''',
        default=None
        )
    parser.add_argument(
        '--gene-families','-f',type=str,
        help='''Provide the path to a file containing the coding sequences for all gene families. This should be a csv matching 
        the output format from Panaroo.''',
        default=None
        )
    parser.add_argument(
        '--rare-genes','-r',type=str,
        help='''Provide the path to a file containing the rare gene families to check. This should be a one-column file.''',
        default=None
        )
    parser.add_argument(
        '--gml','-g',type=str,
        help='''Provide the path to a .gml file from Panaroo that will be used to determine gene adjacency.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output file.''',
        default=None
        )
    args = parser.parse_args()
    check_all_rare_genes(args.locations, args.gene_families, args.rare_genes, args.gml, args.output)

if __name__ == '__main__':
    main()