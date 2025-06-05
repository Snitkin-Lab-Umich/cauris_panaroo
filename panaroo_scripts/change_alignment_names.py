import os
from Bio import SeqIO

def change_names_paralogskip(input_dir,output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            with open(input_dir + fname,'r') as fh_in, open(output_dir + groupname + '.fasta','w') as fh_out:
                seen,towrite = set(),[]
                for record in SeqIO.parse(fh_in,'fasta'):
                    isolatename = record.id.split(';')[0]
                    # THIS WILL SKIP ANY PARALOGS
                    if isolatename not in seen:
                        record.id = isolatename
                        record.name = ''
                        record.description = ''
                        towrite.append(record)
                        seen.add(isolatename)
                SeqIO.write(towrite,fh_out,'fasta')

def change_names(input_dir,output_dir):
    for fname in os.listdir(input_dir):
        if fname.endswith('.aln.fas'):
            groupname = fname.split('.aln.fas')[0]
            with open(input_dir + fname,'r') as fh_in, open(output_dir + groupname + '.fasta','w') as fh_out:
                for line in fh_in:
                    if line.startswith('>'):
                        line = line.split(';')[0] + '\n'
                    _ = fh_out.write(line)

if __name__ == '__main__':
    input_dir = '043025_shortread/panaroo_out_v1/aligned_gene_sequences/'
    output_dir = '043025_shortread/panaroo_out_v1/aligned_gene_sequences_rename/'
    change_names_paralogskip(input_dir,output_dir)
