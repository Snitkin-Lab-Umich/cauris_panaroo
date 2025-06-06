import os
import subprocess
from Bio import SeqIO
import argparse

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
    

def remove_duplicate_entries(gff_in,gff_out_dir,fasta_in,fasta_out_dir,include_ortho):
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
        if include_ortho:
            with open(f'{fasta_in}{isolatename}.proteins.fa','r') as fh_in,open(f'{fasta_out_dir}{isolatename}.proteins.fa','w') as fh_out:
                towrite = []
                for record in SeqIO.parse(fh_in,'fasta'):
                    genename = record.id.split('-T')[0]
                    if genename not in toremove:
                        towrite.append(record)
                SeqIO.write(towrite,fh_out,'fasta')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the name of a directory containing the assemblies and annotation files. The gff files should be in [input]/original_gff,
        and the assemblies should be in [input]/original_nucleotide_fasta''',
        default=None
        )
    parser.add_argument(
        '--ortho','-o',action='store_true',
        help='''Provide the name of a directory containing the assemblies and annotation files. The gff files should be in [input]/original_gff,
        and the assemblies should be in [input]/original_nucleotide_fasta''',
        default=False
        )
    args = parser.parse_args()
    gff_in = f'{args.input}/original_gff/'
    fasta_in = f'{args.input}/original_nucleotide_fasta/'
    fasta_protein_in = f'{args.input}/original_protein_fasta/'
    output_dir1 = f'{args.input}/prokka_gff_withdup/'
    output_dir2 = f'{args.input}/prokka_gff_nodup/'
    output_dir3 = f'{args.input}/protein_fasta_nodup/'
    for p in [output_dir1,output_dir2,output_dir3]:
        if not os.path.exists(p):
            os.mkdir(p)
    convert_all(gff_in,fasta_in,output_dir1)
    remove_duplicate_entries(output_dir1,output_dir2,fasta_protein_in,output_dir3,args.ortho)


if __name__ == '__main__':
    main()