# cauris_panaroo

This is a collection of scripts meant to process a dataset of Candida auris assemblies and generate a pangenome with Panaroo.
If needed, scripts are included to splice introns out of each assembly. This can be done through panaroo_scripts/convert_wrapper_splice_introns.py.
The downstream goal of this script is to run the spliced fasta and gff files through Panaroo's convert_refseq_to_prokka_gff.py script. A copy of this script is included here, originally downloaded from:

https://github.com/gtonkinhill/panaroo
Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. 2020. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21:180.
(See third_party/panaroo/LICENSE)

After processing all assemblies and running Panaroo, it is recommended to run the output through BLAST validation.

