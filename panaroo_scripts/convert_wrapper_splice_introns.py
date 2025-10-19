import os
import subprocess
from Bio import SeqIO
import argparse
import gffutils as gff
import copy

def convert_all(gff_in,fasta_in,output_dir):
    for fname in os.listdir(gff_in):
        fname2 = fname.split('.gff3')[0]
        gff_file = os.path.join(gff_in, fname2 + '.gff3')
        fasta_file = os.path.join(fasta_in, fname2 + '.scaffolds.fa')
        outfile = os.path.join(output_dir, fname2 + '.gff')
        #gff_file = gff_in + fname2 + '.gff3'
        #fasta_file = fasta_in + fname2 + '.scaffolds.fa'
        #outfile = output_dir + fname2 + '.gff'
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
    with open(logfile,'a') as fh_log:
        for fname in os.listdir(gff_input_dir):
            if '.gff' not in fname:
                continue
            isolatename = fname.split('.gff')[0]
            gff_file = os.path.join(gff_input_dir,fname)
            gff_db = gff.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
            removed_intron_count,retained_count = 0,0
            cds_ids = set()
            with open(os.path.join(gff_output_dir,f'{isolatename}.gff3'),'w') as fh_out:
                _ = fh_out.write('##gff-version 3\n')
                #for gene in gff_db.features_of_type('gene'):
                for mrna in gff_db.features_of_type('mRNA'):
                    #gene_id = get_gff_line_ID(gene.attributes['ID'])
                    cds_list = list(gff_db.children(mrna, featuretype='CDS', order_by='start'))
                    if len(cds_list) == 1:
                        # ensure unique ID
                        # note that cds.id will not be the same as cds.attributes['ID'][0]
                        # since merge_strategy='create_unique' was used, cds.id will have a suffix such as _1 or _2 added
                        # cds.attributes['ID'][0] will be the original ID from the gff file, which is the same for all cds features under the same mRNA
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
            _ = fh_log.write(f'Removed {removed_intron_count} genes with multiple CDS features of {total_genes_processed} total genes in {isolatename}\n')


def splice_introns(gff_input_dir,gff_output_dir,fasta_input_dir,fasta_output_dir,logfile):
    # take a gff as the input file
    # check each gene for the number of child CDS features
    # if there's only one, keep that CDS feature
    # if there are exactly two, delete the intron from the fasta file and combine the two CDS features into one
    # if there's more than two, skip that gene entirely, and record the number of genes skipped
    # also ensure that all CDS features have unique IDs
    with open(logfile,'a') as fh_log:
        for fname in os.listdir(gff_input_dir):
            print(f'Processing {fname}')
            if '.gff' not in fname:
                continue
            isolatename = fname.split('.gff')[0]
            # read in gff
            gff_file = os.path.join(gff_input_dir,fname)
            gff_db = gff.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
            # read in corresponding fasta
            fasta_file = os.path.join(fasta_input_dir,f'{isolatename}.scaffolds.fa')
            fasta_records = {record.id:record for record in SeqIO.parse(fasta_file,'fasta')}
            # count spliced, removed, and retained genes
            removed_count,spliced_count,retained_count,isoform_count,framed_count = 0,0,0,0,0
            # keep track of cds ids to ensure uniqueness
            cds_ids = set()
            with open(os.path.join(gff_output_dir,f'{isolatename}.gff3'),'w') as fh_out:
                _ = fh_out.write('##gff-version 3\n')
                # iterate through each scaffold in the gff
                # (the order should not actually matter)
                for scaffold in sorted(list(gff_db.seqids())):
                    # get all genes on this scaffold
                    gene_list = list(gff_db.region(seqid=scaffold,featuretype='gene'))
                    if len(gene_list) == 0:
                        continue
                    # select the longest isoform for each gene
                    mrna_list,isoforms = single_isoform(gene_list,gff_db)
                    isoform_count += isoforms
                    # convert each mRNA to its corresponding CDS features
                    full_cds_list = [list(gff_db.children(mrna, featuretype='CDS', order_by='start')) for mrna in mrna_list]
                    # first, remove any genes that are out of frame
                    cds_list = []
                    for cds_l in full_cds_list:
                        if frame_checker(cds_l):
                            cds_list.append(cds_l)
                        else:
                            framed_count += 1
                    # next, identify all intron regions to be removed from the fasta file, without making any changes
                    initial_intron_regions = find_all_introns_in_cds_list(cds_list)
                    # if there are any cds features in this list with exons that overlap an intron, remove the entire gene
                    cds_list,removed = remove_overlapping_exon_genes(cds_list,initial_intron_regions)
                    # recalculate intron regions based on the remaining entries in cds_list, to account for removed genes with introns
                    intron_regions = find_all_introns_in_cds_list(cds_list)
                    # make a copy of the original cds_list for later verification
                    old_cds_list = copy.deepcopy(cds_list)
                    # then, splice all introns from from the CDS features on this scaffold
                    cds_list,retained,spliced = splice_introns_from_cds_list(cds_list)
                    removed_count += removed
                    retained_count += retained
                    spliced_count += spliced
                    # remove the intron regions from the corresponding fasta record, using the original coordinates
                    if len(intron_regions) > 0:
                        fasta_record = fasta_records[scaffold]
                        old_fasta_record = fasta_record.seq
                        new_fasta_record = fasta_seq_remover(fasta_record,intron_regions)
                        fasta_records[scaffold] = new_fasta_record
                        # verify that the coding sequences without introns have the same exact sequences as before splicing was performed
                        cds_sequence_verifier(old_cds_list,old_fasta_record,cds_list,new_fasta_record.seq,fname)
                    # if all checks were passed, write the modified cds_list to the output gff
                    for cds_l in cds_list:
                        for cds in cds_l:
                            if cds.attributes['ID'][0] in cds_ids:
                                print(f'Error: duplicate CDS ID {cds.attributes["ID"][0]} found in {isolatename}!')
                                quit(1)
                            else:
                                _ = fh_out.write(str(cds) + '\n')
                                cds_ids.add(cds.attributes['ID'][0])
            # write the modified fasta records to a new fasta file
            with open(os.path.join(fasta_output_dir,f'{isolatename}.scaffolds.fa'),'w') as fh_fasta_out:
                SeqIO.write(fasta_records.values(),fh_fasta_out,'fasta')
            total_genes_processed = spliced_count + retained_count + removed_count
            _ = fh_log.write(
                f'Results for {isolatename}:\n'
                f'Resolved {isoform_count} genes with multiple mRNA features. Removed {framed_count} genes with frameshifts. \n' 
                f'Spliced {spliced_count} genes with introns. Retained {retained_count} genes without introns. \n'
                f'Removed {removed_count} genes with exons overlapping introns. {total_genes_processed} total genes processed.\n'
                )


def remove_overlapping_exon_genes(cds_list,intron_regions):
    # take a list of cds features, all of which are on the same scaffold
    # for each entry in the list, compare the cds features to the intron regions
    # if any cds feature overlaps an intron, remove that entry from cds_list
    toremove = []
    remove_count = 0
    for cds_l in cds_list:
        for cds in cds_l:
            for intron_start,intron_end in intron_regions:
                # check if the end of the gene falls within the intron
                if cds.end >= intron_start and cds.end <= intron_end:
                    toremove.append(cds_l)
                    remove_count += 1
                    break
                # check if the start of the gene falls within the intron
                if cds.start >= intron_start and cds.start <= intron_end:
                    toremove.append(cds_l)
                    remove_count += 1
                    break
                # check if the gene completely encompasses the intron
                if cds.start <= intron_start and cds.end >= intron_end:
                    toremove.append(cds_l)
                    remove_count += 1
                    break
    for cds_l in toremove:
        if cds_l in cds_list:
            cds_list.remove(cds_l)
    return(cds_list, remove_count)


def single_isoform(gene_list,gff_db):
    # take a list of gene features from a gffutils database
    # return a list of mRNA features, where the longest isoform for each gene was selected
    mrna_list = []
    isoform_count = 0
    for gene in gene_list:
        mrnas = list(gff_db.children(gene, featuretype='mRNA', order_by='start'))
        if len(mrnas) == 0:
            continue
        elif len(mrnas) == 1:
            mrna_list.append(mrnas[0])
        else:
            # select the longest isoform
            longest_mrna = max(mrnas, key=lambda x: x.end - x.start)
            mrna_list.append(longest_mrna)
            isoform_count += 1
    return(mrna_list,isoform_count)

def cds_sequence_verifier(old_cds_list,old_fasta_record,new_cds_list,new_fasta_record,fname):
    # for each cds in old_cds_list, if it does not have an intron, get its nucleotide sequence from old_fasta_record
    # then, verify that the corresponding cds in new_cds_list has the same nucleotide sequence from new_fasta_record
    for ind in range(len(old_cds_list)):
        old_cds = old_cds_list[ind]
        if len(old_cds) != 1:
            continue
        old_cds = old_cds[0]
        # the old and new cds_lists should be in the same order
        new_cds = new_cds_list[ind][0]
        old_cds_seq = old_fasta_record[(old_cds.start - 1):(old_cds.end)]  # convert to 0-based
        new_cds_seq = new_fasta_record[(new_cds.start - 1):(new_cds.end)]  # convert to 0-based
        if old_cds_seq != new_cds_seq:
            print(f'New CDS sequence does not match old CDS sequence for {old_cds.id} in {fname}!')
            print(old_cds.id, old_cds.start, old_cds.end)
            print(old_cds_seq)
            print(new_cds.id, new_cds.start, new_cds.end)
            print(new_cds_seq)
            quit(1)

def fasta_seq_remover(fasta_record,intron_regions):
    # fasta_record is a Biopython SeqRecord object for a single scaffold
    # intron_regions is a list of tuples, where each tuple is (start,end) for a region to remove
    # return a modified fasta_record with the intron regions removed
    # first, sort the intron regions by start coordinate, in descending order
    intron_regions = sorted(intron_regions, key=lambda x: x[0], reverse=True)
    # next, ensure no intron regions overlap, taking into account that these are in descending order
    # compare the start of the downstream intron [i] to the end of the upstream intron [i+1]
    for i in range(len(intron_regions)-1):
        if intron_regions[i][0] <= intron_regions[i+1][1]:
            print(intron_regions[i], intron_regions[i+1])
            print('Overlapping introns detected!')
            quit(1)
    # finally, remove each intron region from the fasta record
    for intron_interval in intron_regions:
        intron_start = intron_interval[0] - 1  # convert to 0-based
        intron_end = intron_interval[1]  # end is exclusive in python slicing
        fasta_record.seq = fasta_record.seq[:intron_start] + fasta_record.seq[intron_end:]
    return(fasta_record)

def find_all_introns_in_cds_list(cds_list):
    # take a list of genes that have already been converted to their child cds features
    # all cds features should be on the same scaffold/seqid
    # return a list of intron regions to be removed from the fasta file
    # do not modify the cds_list
    intron_regions = []
    for cds_l in cds_list:
        if len(cds_l) > 1:
            # for each pair of adjacent cds features, find the intron between them and record it
            for i in range(len(cds_l)-1):
                upstream_cds = cds_l[i]
                downstream_cds = cds_l[i+1]
                intron_start,intron_end,intron_length = intron_finder(upstream_cds,downstream_cds)
                intron_regions.append((intron_start,intron_end))
    return intron_regions

def splice_introns_from_cds_list(cds_list):
    # take a list of genes that have already been converted to their child cds features
    # all cds features should be on the same scaffold/seqid
    # return a modified cds_list, where introns have been spliced out of genes with exactly two cds features
    # also modify the coordinates of all downstream features whenever an intron is removed
    # finally, return a list of intron regions to be removed from the fasta file
    retained, spliced, framed = 0, 0, 0
    for cds_l_index in range(len(cds_list)):
        cds_l = cds_list[cds_l_index]
        if len(cds_l) == 1:
            retained += 1
        elif len(cds_l) > 1:
            spliced += 1
            # as long as at least one intron remains:
            while len(cds_l) > 1:
                # take the first and second coding sequences in the list
                # the list should be sorted by start position, and consist entirely of cds features from the same mRNA, so these should not overlap
                upstream_cds = cds_l[0]
                downstream_cds = cds_l[1]
                # locate the intron between these two coding sequences
                intron_start,intron_end,intron_length = intron_finder(upstream_cds,downstream_cds)
                # modify cds_list to combine these two introns
                cds_list = intron_remover(cds_list,cds_l_index,intron_end,intron_length)
                # update cds_l to reflect the modified cds_list
                cds_l = cds_list[cds_l_index]
            # after all introns have been removed, update the frame attribute to '0'
            cds_l[0].frame = '0'
    return(cds_list,retained,spliced)

def frame_checker(cds_l):
    # take a list of cds features from a single gene, sorted by start position
    # verify that this is not a frame shift with two checks:
    # the total length of all cds features should be a multiple of 3
    # the first or last cds feature should have a frame attribute of '0'
    total_length = sum([cds.end - cds.start + 1 for cds in cds_l])
    if total_length % 3 != 0:
        return(False)
    elif cds_l[0].frame != '0' and cds_l[-1].frame != '0':
        return(False)
    else:
        return(True)

def intron_finder(upstream_cds,downstream_cds):
    # given two cds annotations, return the coordinates of the intron between them
    # first, verify that the two cds features are in the expected order and format
    if not cds_checker(upstream_cds,downstream_cds):
        print(f'Error: unexpected CDS features for {upstream_cds.id} and {downstream_cds.id}')
        quit(1)
    intron_start = upstream_cds.end + 1
    intron_end = downstream_cds.start - 1
    intron_length = intron_end - intron_start + 1
    return(intron_start,intron_end,intron_length)

def cds_checker(upstream_cds,downstream_cds):
    # verify that two cds features are in the expected order and format
    checks = [upstream_cds.start > upstream_cds.end]
    checks.append(downstream_cds.start > downstream_cds.end)
    checks.append(upstream_cds.end > downstream_cds.start)
    checks.append(upstream_cds.start > downstream_cds.start)
    checks.append(upstream_cds.strand != downstream_cds.strand)
    if any(checks):
        return(False)
    else:
        return(True)

def intron_remover(cds_list,cds_l_index,intron_end,intron_length):
    # cds_list is a list of genes converted to their child cds features
    # cds_l_index is the index of the cds_list containing the intron to be removed
    # (this cds_list contains exactly two cds features)
    # first, adjust all coordinates in cds_list that are downstream of the intron
    for cds_l in cds_list:
        for cds in cds_l:
            if cds.start > intron_end:
                cds.start -= intron_length
                cds.end -= intron_length
    # make sure the two cds features are now directly adjacent
    upstream_cds = cds_list[cds_l_index][0]
    downstream_cds = cds_list[cds_l_index][1]
    if upstream_cds.end + 1 != downstream_cds.start:
        print(f'Error: intron removal did not result in adjacent CDS features for {upstream_cds.id} and {downstream_cds.id}')
        quit(1)
    # combine the two cds features into one
    new_cds = upstream_cds
    new_cds.end = downstream_cds.end
    # replace the upstream cds with the new cds, and remove the downstream cds
    cds_list[cds_l_index][0] = new_cds
    cds_list[cds_l_index].remove(downstream_cds)
    return(cds_list)

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

def copy_data(isolate,annotate_results_dir,gff_input,nucl_fasta_input,prot_fasta_input):
    if not os.path.isdir(annotate_results_dir):
        print(f'Could not locate annotation results for {isolate} at {annotate_results_dir}')
        return
    for suff,dest in [['.gff3',gff_input],['.scaffolds.fa',nucl_fasta_input],['.proteins.fa',prot_fasta_input]]:
        source_file = os.path.join(annotate_results_dir,f'{isolate}{suff}')
        if not os.path.isfile(source_file):
            print(f'Could not locate expected file {source_file}')
            return
        subprocess.run(['cp',source_file,dest])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--assemblies','-a',type=str,
        help='''Provide the name of a directory containing funannotate outputs for all isolates. This should be in the same format as 
        cauris-data-flow's funannotate directory.''',
        default=None
        )
    parser.add_argument(
        '--additional-assemblies','-aa',type=str,
        help='''Provide the name of a directory containing funannotate outputs for all isolates. This should be in the same format as 
        cauris-data-flow's funannotate directory. This option is intended for adding hybrid assemblies, and will overwrite the files provided 
        in --assemblies if they share the same name.''',
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
    parser.add_argument(
        '--introns','-i',choices=['remove','splice'],type=str,
        help='''Specify whether to remove genes with introns or attempt to splice out introns. Default is removal. Genes with multiple
        introns will be removed regardless of the option chosen.''',
        default='remove'
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
            isolate_annotate_results = os.path.join(args.assemblies,isolate,'annotate_results')
            copy_data(isolate,isolate_annotate_results,gff_input,nucl_fasta_input,prot_fasta_input)
        if args.additional_assemblies is not None:
            if not os.path.isdir(args.additional_assemblies):
                print(f'Could not locate input directory at {args.additional_assemblies}')
                quit(1)
            for isolate2 in os.listdir(args.additional_assemblies):
                isolate_additional_annotate_results = os.path.join(args.additional_assemblies,isolate2,'annotate_results')
                copy_data(isolate2,isolate_additional_annotate_results,gff_input,nucl_fasta_input,prot_fasta_input)
        print(f'Finished copying files from {args.assemblies} to output directory {args.output}')
    gff_introns_removed = os.path.join(args.output,f'original_gff_{args.introns}Introns')
    gff_final = os.path.join(args.output,f'prokka_gff_{args.introns}Introns')
    logfile = os.path.join(log_dir,f'{args.introns}Introns_log.txt')
    for p in [gff_introns_removed,gff_final,log_dir]:
        if not os.path.exists(p):
            os.mkdir(p)
    # write a single line to the log file to initialize
    with open(logfile,'w') as fh_log:
        _ = fh_log.write(f'Intron and duplicate removal log - {args.introns} mode\n')
    if args.introns == 'remove':
        remove_intron_genes(gff_input,gff_introns_removed,logfile)
    elif args.introns == 'splice':
        fasta_introns_removed = os.path.join(args.output,f'nucleotide_fasta_{args.introns}Introns')
        if not os.path.exists(fasta_introns_removed):
            os.mkdir(fasta_introns_removed)
        splice_introns(gff_input,gff_introns_removed,nucl_fasta_input,fasta_introns_removed,logfile)
        nucl_fasta_input = fasta_introns_removed
    convert_all(gff_introns_removed,nucl_fasta_input,gff_final)


if __name__ == '__main__':
    # gff_input = 'test/original_gff/'
    # gff_introns_removed = 'test/spliced_gff/'
    # nucl_fasta_input = 'test/original_fasta/'
    # fasta_introns_removed = 'test/spliced_fasta/'
    # logfile = 'test/splice_log.txt'
    # splice_introns(gff_input,gff_introns_removed,nucl_fasta_input,fasta_introns_removed,logfile)
    main()