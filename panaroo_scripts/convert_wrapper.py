import os
import subprocess
from Bio import SeqIO
import argparse
import gffutils as gff

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


def remove_intron_genes(gff_input_dir,gff_output_dir,logfile):
    # take a gff as the input file
    # check each gene for the number of child CDS features
    # if there's only one, keep that CDS feature
    # if there's more than one, skip that gene entirely, and record the number of genes skipped
    # also ensure that all CDS features have unique IDs
    with open(logfile,'w') as fh_log:
        for fname in os.listdir(gff_input_dir):
            if '.gff' not in fname:
                continue
            isolatename = fname.split('.gff')[0]
            gff_db = gff.create_db(fh_in, dbfn=':memory:', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
            removed_intron_count,retained_count = 0,0
            cds_ids = set()
            with open(f'{gff_output_dir}{isolatename}.gff3','w') as fh_out:
                _ = fh_out.write('##gff-version 3\n')
                for gene in gff_db.features_of_type('gene'):
                    #gene_id = get_gff_line_ID(gene.attributes['ID'])
                    cds_list = list(gff_db.children(gene, featuretype='CDS', order_by='start'))
                    if len(cds_list) == 1:
                        # ensure unique ID
                        cds = cds_list[0]
                        if cds.attributes['ID'][0] in cds_ids:
                            print(f'Error: duplicate CDS ID {cds.attributes["ID"][0]} found in {isolatename}!')
                            quit(1)
                        else:
                            _ = fh_out.write(str(cds_list[0]) + '\n')
                            cds_ids.add(cds.attributes['ID'][0])
                            retained_count += 1
                    else:
                        removed_intron_count += 1
            total_genes_processed = removed_intron_count + retained_count
            _ = fh_log.write(f'Removed {removed_intron_count} genes with multiple CDS features of {total_cds_processed} total genes in {isolatename}')


def splice_introns(gff_input_dir,gff_output_dir):
    # take a gff as the input file
    # look at each gene and its child CDS features
    # if that gene has multiple child CDS features, combine them into a single CDS feature that spans from the start of the first to the end of the last
    # if there is only one CDS feature, just write it out as is
    # write out a new gff file that only has one CDS feature per gene
    for fname in os.listdir(gff_input_dir):
        isolatename = fname.split('.gff')[0]
        gff_db = gff.create_db(fh_in, dbfn=':memory:', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
        genes = {}
        for gene in gff_db.features_of_type('gene'):
            gene_id = get_gff_line_ID(gene.attributes['ID'])
            genes[gene_id] = []
            for cds in gff_db.children(gene, featuretype='CDS', order_by='start'):
                genes[gene_id].append(cds)
        with open(f'{gff_output_dir}{fname}','w') as fh_out, open(log,'w') as fh_log:
            _ = fh_out.write('##gff-version 3\n')
            for gene_id in genes:
                if len(genes[gene_id]) == 1:
                    _ = fh_out.write(str(cds) + '\n')
                if len(genes[gene_id]) > 1:
                    _ = fh_log.write(f'{len(genes[gene_id])} CDS features for {gene_id} were collapsed\n')
                    # combine the CDS features
                    chrom = genes[gene_id][0].chrom
                    source = genes[gene_id][0].source
                    featuretype = 'CDS'
                    start = min([cds.start for cds in genes[gene_id]])
                    end = max([cds.end for cds in genes[gene_id]])
                    score = genes[gene_id][0].score
                    strand = genes[gene_id][0].strand
                    phase = genes[gene_id][0].frame
                    attributes = genes[gene_id][0].attributes
                    new_cds = gff.Feature(chrom, source, featuretype, start, end, score, strand, phase, attributes)
                    _ = fh_out.write(str(new_cds) + '\n')
                    _ = fh_log.write(f'Collapsed {len(genes[gene_id])} CDS features for gene:{gene_id}\n')



def remove_duplicate_entries(gff_in,gff_out_dir,fasta_in,fasta_out_dir,include_ortho,logfile):
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
        '--assemblies','-a',type=str,
        help='''Provide the name of a directory containing funannotate outputs for all isolates. This should be in the same format as 
        cauris-data-flow's funannotate directory.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the name of a directory where the output files should go. A large number of files will be created here.''',
        default=None
        )
    parser.add_argument(
        '--skip-copy','-sc',action='store_true',
        help='''Set this flag to skip copying the original funannotate files to the output directory. The script will assume that original files are already in
        the output directory. Processed gff files will still be created.''',
        default=False
        )
    args = parser.parse_args()
    # the expected format of the input directory is:
    # funannotate/[isolate_name]/annotate_results/
    # this directory should contain:
    # [isolate_name].scaffolds.fa
    # [isolate_name].proteins.fa
    # [isolate_name].gff3
    # all of these need to be copied to their respective subdirectories in the output directory
    gff_input = os.path.join(args.output,'original_gff')
    nucl_fasta_input = os.path.join(args.output,'original_nucleotide_fasta')
    prot_fasta_input = os.path.join(args.output,'original_protein_fasta')
    log_dir = os.path.join(args.output,'logs')
    if not args.skip_copy:
        if not os.path.isdir(args.assemblies):
            print(f'Could not locate input directory at {args.assemblies}')
            quit(1)
        if not os.path.exists(args.output):
            os.mkdir(args.output)
        for p in [gff_input,nucl_fasta_input,prot_fasta_input]:
            if not os.path.exists(p):
                os.mkdir(p)
        for isolate in os.listdir(args.assemblies):
            isolate_path = os.path.join(args.assemblies,isolate,'annotate_results')
            if not os.path.isdir(isolate_path):
                print(f'Could not locate annotation results for isolate {isolate} at {isolate_path}')
                continue
            for suff,dest in [['.gff3',gff_input],['.scaffolds.fa',nucl_fasta_input],['.proteins.fa',prot_fasta_input]]:
                source_file = os.path.join(isolate_path,f'{isolate}{suff}')
                if not os.path.isfile(source_file):
                    print(f'Could not locate expected file {source_file}')
                    continue
                subprocess.run(['cp',source_file,dest])
        print(f'Finished copying files from {args.assemblies} to output directory {args.output}')
    gff_introns_removed = os.path.join(args.output,'original_gff_NoIntrons')
    gff_final = os.path.join(args.output,'prokka_gff')
    logfile = os.path.join(log_dir,'intron_removal_log.txt')
    for p in [gff_introns_removed,gff_final,log_dir]:
        if not os.path.exists(p):
            os.mkdir(p)
    remove_intron_genes(gff_input,gff_introns_removed,logfile)
    convert_all(gff_introns_removed,nucl_fasta_input,gff_final)


if __name__ == '__main__':
    main()