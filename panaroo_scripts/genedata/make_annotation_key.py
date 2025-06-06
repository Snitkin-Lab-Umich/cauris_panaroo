
import os
import gffutils as gff
from Bio import SeqIO

def table_to_dict(file_path):
    outdict = {}
    title_read = False
    with open(file_path,'r') as fh:
        for line in fh:
            if not title_read:
                keylist = line.strip().split('\t')
                title_read = True
                continue
            llist = line.strip().split('\t')
            ddict = {}
            for i in range(len(keylist)):
                ddict[keylist[i]] = llist[i]
            outdict[ddict['mRNAID']] = ddict
    return(outdict)

def make_key_from_ref(reftable,querytable,query_isolate,outfile):
    ref_dict = table_to_dict(reftable)
    q_dict = table_to_dict(querytable)
    # this makes two dictionaries with the mRNA ID as a key
    # all values are dictionaries, and these dictionaries should have the same keys
    with open(outfile,'a') as fh:
        for q in q_dict:
            q_entry = q_dict[q]
            #print(q_entry)
            for r in ref_dict:
                ref_entry = ref_dict[r]
                #print(ref_entry)
                # for now, only look at well-annotated genes
                tests = [q_entry['Name'] == ref_entry['Name']]
                tests.append(q_entry['Name'] != 'NA')
                tests.append(q_entry['Product'] == ref_entry['Product'])
                tests.append(q_entry['Product'] != 'hypothetical protein')
                if all(tests):
                    tests2 = 0
                    for k in ['BUSCO','COG','EggNog','IPR','PFAM']:
                        if q_entry[k] == ref_entry[k]:
                            tests2+=1
                    for c in ['CDSLength','GeneLength']:
                        if q_entry[c] != 'NA' and ref_entry[c] != 'NA':
                            if abs(int(q_entry[c]) - int(ref_entry[c])) < 300:
                                tests2+=1
                    if tests2 > 5:
                        line = '\t'.join([query_isolate + '_' + q_entry['mRNAID'],query_isolate + '_' + q_entry['Name']]) + '\n'
                        _ = fh.write(line)

def main(ref_table,query_dir,outfile):
    querylist = [x for x in os.listdir(query_dir) if x.endswith('.tsv')]
    # copy header line from ref table
    with open(outfile,'w') as fh_out:
        _  = fh_out.write('original_mRNAID\tnew_ID\n')
    for q in querylist:
        query_isolate = q.split('.tsv')[0]
        query_table = query_dir + q
        if query_table != ref_table:
            make_key_from_ref(ref_table,query_table,query_isolate,outfile)

if __name__ == '__main__':
    ref_table = 'panaroo_scripts/genedata/original_gff_annotation_data/UM_Caur_1.tsv'
    query_dir = 'panaroo_scripts/genedata/original_gff_annotation_data/'
    outfile = 'panaroo_scripts/genedata/conversion_table_v1.tsv'
    main(ref_table,query_dir,outfile)


# def clean_gff_string(gff_string):
#     splitlines = gff_string.splitlines()
#     lines_to_delete = []
#     for index in range(len(splitlines)):
#         if '##sequence-region' in splitlines[index]:
#             lines_to_delete.append(index)
#     for index in sorted(lines_to_delete, reverse=True):
#         del splitlines[index]
#     cleaned_gff = "\n".join(splitlines)
#     return cleaned_gff


# def convert(gfffile, outputfile, fastafile, is_ignore_overlapping):

#     #Split file and parse
#     with open(gfffile, 'r') as infile:
#         lines = infile.read().replace(',','')

#     if fastafile is None:
#         split = lines.split('##FASTA')
#         if len(split) != 2:
#             print("Problem reading GFF3 file: ", gfffile)
#             raise RuntimeError("Error reading GFF3 input!")
#     else:
#         with open(fastafile, 'r') as infile:
#             fasta_lines = infile.read()
#         split = [lines, fasta_lines]

#     with StringIO(split[1]) as temp_fasta:
#         sequences = list(SeqIO.parse(temp_fasta, 'fasta'))

#     for seq in sequences:
#         seq.description = ""

#     parsed_gff = gff.create_db(clean_gff_string(split[0]),
#                                dbfn=":memory:",
#                                force=True,
#                                keep_order=False,
#                                merge_strategy="create_unique",
#                                sort_attribute_values=True,
#                                from_string=True)

#     with open(outputfile, 'w') as outfile:
#         # write gff part
#         outfile.write("##gff-version 3\n")
#         for seq in sequences:
#             outfile.write(
#                 " ".join(["##sequence-region", seq.id, "1",
#                           str(len(seq.seq))]) + "\n")

#         prev_chrom = ""
#         prev_end = -1
#         ids = set()
#         seen = set()
#         seq_order = []
#         for entry in parsed_gff.all_features(featuretype=(),
#                                              order_by=('seqid', 'start')):
#             entry.chrom = entry.chrom.split()[0]
#             # skip non CDS
#             if "CDS" not in entry.featuretype: continue
#             # skip overlapping CDS if option is set
#             if entry.chrom == prev_chrom and entry.start < prev_end and is_ignore_overlapping:
#                 continue
#             # skip CDS that dont appear to be complete or have a premature stop codon

#             premature_stop = False
#             for sequence_index in range(len(sequences)):
#                 scaffold_id = sequences[sequence_index].id
#                 if scaffold_id == entry.seqid:
#                     gene_sequence = sequences[sequence_index].seq[(
#                         entry.start - 1):entry.stop]
#                     if (len(gene_sequence) % 3 > 0) or (len(gene_sequence) <
#                                                         34):
#                         premature_stop = True
#                         break
#                     if entry.strand == "-":
#                         gene_sequence = gene_sequence.reverse_complement()
#                     if "*" in str(gene_sequence.translate())[:-1]:
#                         premature_stop = True
#                         break
#             if premature_stop: continue

#             c = 1
#             while entry.attributes['ID'][0] in ids:
#                 entry.attributes['ID'][0] += "." + str(c)
#                 c += 1
#             ids.add(entry.attributes['ID'][0])
#             prev_chrom = entry.chrom
#             prev_end = entry.end
#             if entry.chrom not in seen:
#                 seq_order.append(entry.chrom)
#                 seen.add(entry.chrom)
#             print(entry, file=outfile)

#         # write fasta part
#         outfile.write("##FASTA\n")
#         sequences = [
#             seq for x in seq_order for seq in sequences if seq.id == x
#         ]
#         if len(sequences) != len(seen):
#             raise RuntimeError("Mismatch between fasta and GFF!")
#         SeqIO.write(sequences, outfile, "fasta")

#     return


# def main():

#     parser = argparse.ArgumentParser(
#         description='Converts refseq GFF3 to prokka format.')

#     parser.add_argument('-g',
#                         '--gff',
#                         dest='gff',
#                         type=str,
#                         required=True,
#                         help='input gff file name')

#     parser.add_argument(
#         '-f',
#         '--fasta',
#         dest='fasta',
#         type=str,
#         default=None,
#         help='input fasta file name (if separate from the GFF)')

#     parser.add_argument('-o',
#                         '--out',
#                         dest='out',
#                         type=str,
#                         required=True,
#                         help='output file name')
#     parser.add_argument(
#         '--is_ignore_overlapping',
#         action="store_true",
#         help="set to ignore CDS that overlap (that's common in bacteria)")
#     args = parser.parse_args()

#     convert(args.gff, args.out, args.fasta, args.is_ignore_overlapping)

#     return


# if __name__ == '__main__':
#     main()
