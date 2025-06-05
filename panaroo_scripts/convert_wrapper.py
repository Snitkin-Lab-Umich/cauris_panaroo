import os
import subprocess
from Bio import SeqIO

def convert_all(gff_in,fasta_in,output_dir):
    for fname in os.listdir(gff_in):
        fname2 = fname.split('.gff3')[0]
        gff_file = gff_in + fname2 + '.gff3'
        fasta_file = fasta_in + fname2 + '.scaffolds.fa'
        outfile = output_dir + fname2 + '.gff'
        command = ['python3.12','panaroo_scripts/convert_refseq_to_prokka_gff.py','--gff',gff_file,'--fasta',fasta_file,'--out',outfile]
        subprocess.run(command)

def get_gff_line_ID(linelist):
    ids = linelist[8]
    if 'ID=' not in ids:
        print(f'error with: {ids}')
        quit(1)
    cdsname = ids.split(';')[0].split('ID=')[1]
    # isolate all text between ID= and ;
    # this should always be in the format FUN_000096-T1.cds
    cdsname = cdsname.split('-T')[0]
    return(cdsname)
    

def remove_duplicate_entries(gff_in,gff_out_dir,fasta_in,fasta_out_dir):
    for fname in os.listdir(gff_in):
        isolatename = fname.split('.gff')[0]
        with open(f'{gff_in}{fname}','r') as fh_in,open(f'{gff_out_dir}{fname}','w') as fh_out:
            n=set()
            toremove=set()
            for line in fh_in:
                if line.startswith('scaffold_'):
                    # all annotation lines should start with this
                    line2=line.strip().split('\t')
                    location=line2[0]
                    start=min(int(line2[3]),int(line2[4]))
                    end=max(int(line2[3]),int(line2[4]))
                    genename = get_gff_line_ID(line2)
                    d=(location,start,end)
                    if d in n:
                        # any annotation with the exact same coordinates is added to a list for removal from orthofinder's fastas as well
                        toremove.add(genename)
                        continue
                    n.add(d)
                _ = fh_out.write(line)
        # with this info, remove the same entries from the orthofinder fastas
        with open(f'{fasta_in}{isolatename}.proteins.fa','r') as fh_in,open(f'{fasta_out_dir}{isolatename}.proteins.fa','w') as fh_out:
            towrite = []
            for record in SeqIO.parse(fh_in,'fasta'):
                genename = record.id.split('-T')[0]
                if genename not in toremove:
                    towrite.append(record)
            SeqIO.write(towrite,fh_out,'fasta')
                



if __name__ == '__main__':
    gff_in = '043025_shortread/original_gff/'
    fasta_in = '043025_shortread/original_nucleotide_fasta/'
    fasta_protein_in = '043025_shortread/original_protein_fasta/'
    output_dir1 = '043025_shortread/prokka_gff_withdup/'
    output_dir2 = '043025_shortread/prokka_gff_nodup/'
    output_dir3 = '043025_shortread/protein_fasta_nodup/'
    for p in [output_dir1,output_dir2,output_dir3]:
        if not os.path.exists(p):
            os.mkdir(p)
    convert_all(gff_in,fasta_in,output_dir1)
    remove_duplicate_entries(output_dir1,output_dir2,fasta_protein_in,output_dir3)
