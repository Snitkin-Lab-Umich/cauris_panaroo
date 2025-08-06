
import os
import gffutils as gff
from Bio import SeqIO

def read_fasta(fasta_file):
    # convert a fasta file into a python dictionary
    d = {}
    with open(fasta_file,'r') as fh:
        for record in SeqIO.parse(fh,'fasta'):
            d[record.id] = str(record.seq)
    return(d)

def searchlist(term,datalist):
    termlist = [x.replace(term,'') for x in datalist if term in x]
    if any(termlist):
        return(','.join(sorted(termlist)))
    else:
        return(None)


def get_attribute_data(feature,geneID_to_name):
    attrib_data = {
        'GeneID':None,'mRNAID':None,'CDSID':None,'Name':None,'Product':None,'IPR':None,'BUSCO':None,
        'PFAM':None,'COG':None,'EggNog':None,'GeneLength':None,'CDSLength':None,'CDS_GC':None}
    # add the mRNA ID and gene ID (the gene ID should be the parent)
    mrna_id = ','.join(feature.attributes['ID'])
    gene_id = ','.join(feature.attributes['Parent'])
    attrib_data['mRNAID'] = mrna_id
    attrib_data['GeneID'] = gene_id
    # check dict for name (all mRNA features should have one parent)
    attrib_data['Name'] = geneID_to_name.get(gene_id)
    # get length of feature - note that this will be the length of the entire annotation, not just the CDS
    attrib_data['GeneLength'] = abs(feature.end - feature.start)
    # get product, if present
    if 'product' in feature.attributes.keys():
        attrib_data['Product'] = ','.join(feature.attributes['product'])
    # parse Dbxref (InterProScan and PFAM)
    if 'Dbxref' in feature.attributes.keys():
        dbx_data = feature.attributes['Dbxref']
        attrib_data['IPR'] = searchlist('InterPro:',dbx_data)
        attrib_data['PFAM'] = searchlist('PFAM:',dbx_data)
    if 'note' in feature.attributes.keys():
        note_data = feature.attributes['note']
        attrib_data['EggNog'] = searchlist('EggNog:',note_data)
        attrib_data['COG'] = searchlist('COG:',note_data)
        attrib_data['BUSCO'] = searchlist('BUSCO:',note_data)
    return(attrib_data)

def get_gc_content(fasta_seq,seq_id,seq_start,seq_end):
    # seq start and end may not be in order
    abs_seq_start = (min([seq_start,seq_end]))
    abs_seq_end = (max([seq_start,seq_end]))
    gene_seq = fasta_seq[seq_id][(abs_seq_start-1):(abs_seq_end)]
    gc_count = sum([x in ['G','C'] for x in gene_seq])
    gc_percent = round((gc_count/len(gene_seq)) * 100,1)
    return(gc_percent)

def read_gff(gff_file,fasta_file):
    all_attrib_data = {}
    fasta_seq = read_fasta(fasta_file)
    parsed_gff = gff.create_db(gff_file,
                               dbfn=":memory:",
                               force=True,
                               keep_order=False,
                               merge_strategy="create_unique",
                               sort_attribute_values=True,
                               from_string=False)
    # for each mRNA, extract the annotation info
    # note that this assumes every annotation has a parent, in this order: gene -> mRNA -> CDS
    geneID_to_name = {}
    saved_cds,skipped_cds = 0,0
    for feature in parsed_gff.all_features(featuretype=(),order_by=('seqid', 'start')):
        if 'gene' in feature.featuretype:
            # get the name of this gene and add it to the dict
            gene_name = feature.attributes.get('Name')
            if gene_name is not None:
                gene_name = gene_name[0]
            geneID_to_name[feature.attributes['ID'][0]] = gene_name
        if 'mRNA' in feature.featuretype:
            # get the attribute data of this mRNA, add it to the dictionary
            # use the mRNA name as the key
            attrib_data = get_attribute_data(feature,geneID_to_name)
            all_attrib_data[attrib_data['mRNAID']] = attrib_data            
        if 'CDS' in feature.featuretype:
            cds_id = ','.join(feature.attributes['ID'])
            parent_mrna_id = ','.join(feature.attributes['Parent'])
            if parent_mrna_id in all_attrib_data:
                if all_attrib_data[parent_mrna_id]['CDSID'] is None:
                    # this will skip additional CDS annotations for the same gene/mRNA
                    all_attrib_data[parent_mrna_id]['CDSID'] = cds_id
                    all_attrib_data[parent_mrna_id]['CDSLength'] = abs(feature.end - feature.start)
                    all_attrib_data[parent_mrna_id]['CDS_GC'] = get_gc_content(fasta_seq,feature.seqid,feature.start,feature.end)
                    saved_cds+=1
                else:
                    skipped_cds+=1
    print(f'Skipped {round(skipped_cds/(skipped_cds+saved_cds)*100,1)}% of CDS annotations for {gff_file}')
    return(all_attrib_data)

def write_attrib_data(all_attrib_data,outfile):
    header_write = False
    with open(outfile,'w') as fh:
        for mrna_id in all_attrib_data:
            if not header_write:
                col_order = sorted(list(all_attrib_data[mrna_id].keys()))
                _ = fh.write('\t'.join(col_order) + '\n')
                header_write = True
            data = [all_attrib_data[mrna_id][x] for x in col_order]
            # this puts the contents of all_attrib_data[mrna_id] in the expected order
            data2 = [str(x) if x is not None else 'NA' for x in data]
            # this replaces any None elements with a string instead, and converts numerics into strings
            _ = fh.write('\t'.join(data2) + '\n')

def main(gff_dir,fasta_dir,out_dir,gff_suffix,fasta_suffix):
    namelist = [x.split('.gff')[0] for x in os.listdir(gff_dir) if x.endswith(('.gff','.gff3'))]
    for n in namelist:
        data_dict = read_gff(
            gff_file = gff_dir + n + gff_suffix, 
            fasta_file = fasta_dir + n + fasta_suffix)
        write_attrib_data(all_attrib_data = data_dict, outfile = out_dir + n + '.tsv')




if __name__ == '__main__':
    gff_dir = '040825_shortread/original_gff/'
    gff_suffix = '.gff3'
    fasta_dir = '040825_shortread/original_nucleotide_fasta/'
    fasta_suffix = '.scaffolds.fa'
    out_dir = 'panaroo_scripts/genedata/original_gff_annotation_data/'
    main(gff_dir,fasta_dir,out_dir,gff_suffix,fasta_suffix)

