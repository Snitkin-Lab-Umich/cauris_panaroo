import os
import subprocess

def convert_all(gff_in,fasta_in,output_dir):
    for fname in os.listdir(gff_in):
        fname2 = fname.split('.gff3')[0]
        gff_file = gff_in + fname2 + '.gff3'
        fasta_file = fasta_in + fname2 + '.scaffolds.fa'
        outfile = output_dir + fname2 + '.gff'
        command = ['python3.12','panaroo_scripts/convert_refseq_to_prokka_gff.py','--gff',gff_file,'--fasta',fasta_file,'--out',outfile]
        subprocess.run(command)

def remove_duplicate_entries(gff_in,out_dir):
    for fname in os.listdir(gff_in):
        with open(f'{gff_in}{fname}','r') as fh_in,open(f'{out_dir}{fname}','w') as fh_out:
            n=set()
            for line in fh_in:
                if line.startswith('scaffold_'):
                    # all annotation lines should start with this
                    line2=line.strip().split('\t')
                    location=line2[0]
                    start=min(int(line2[3]),int(line2[4]))
                    end=max(int(line2[3]),int(line2[4]))
                    d=(location,start,end)
                    if d in n:
                        continue
                    n.add(d)
                _ = fh_out.write(line)

if __name__ == '__main__':
    gff_in = '043025_shortread/original_gff/'
    fasta_in = '043025_shortread/original_nucleotide_fasta/'
    output_dir1 = '043025_shortread/prokka_gff_withdup/'
    output_dir2 = '043025_shortread/prokka_gff_nodup/'
    for p in [output_dir1,output_dir2]:
        if not os.path.exists(p):
            os.mkdir(p)
    convert_all(gff_in,fasta_in,output_dir1)
    remove_duplicate_entries(output_dir1,output_dir2)
