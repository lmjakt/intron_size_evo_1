# Scripts and data files for 'Intron size minimisation in teleosts'

This analysis was carried out using a combination of Perl, R and R extensions
written in C. To complicate matters the analysis was carried out on two
separate computers due to operational reasons. An initial R-session was later
split into multiple more clearly defined sessions, but using data objects
exported from the initial session in order to avoid redoing the major parts of
the analysis. Nevertheless, there is an unfortunate amount of overlap between
the initial file and the subsequent targetted analyses.

The files in this git repository come from data extraction and analyses which
were carried out on computer 1. Computer 1 contained the Ensembl databases but
did not have sufficient memory or compute resources for the stages of analysis
depending on large scale sequence alignments.

**[note to self: this is missing the first basic statistics which used all
intron sizes from all genes]**

1. Installation of the majority of Ensembl vertebrate databases using:
     a. `core_dbs/install_db.h`
     b. `core_dbs/install_ensembl_core_98.sh`
2. Extraction of summary statistics and sequences from all canonical
   transcripts using `compile_all_seqs.pl` called from
   `compile_all_seqs.sh`. This produced
   `family_members/all_genes/exon_intron_stats.csv` which was used in
   downstream analyses (see files and process for computer 2).
2. Definition of a suitable set of gene orthologues from which to derive
   intron orthologues in:  
   `family_members/family_members.R`  
   Exported to `family_members/vertebrate_family_members_1_130.txt`
3. Extraction of intron and exon sequences from the full genome
   sequences using a Perl script:  
   `family_members/compile_seqs.pl`  
   Called from `family_members/compile_seqs.sh`.  
   Sequence properties exported to :
   `family_members/orthologue_transcripts/exon_intron_stats.csv`.
4. Preliminary data exploration and export (`R_172_genomes/R_i72_genomes.R`)
   of data structures directly obtained from the core Ensembl databases: 
     a. `R_172_genomes/genome_sizes.txt`
	 b. `R_172_genomes/species_lineages.rds` NCBI taxonomy obtained from
     Ensembl.
	 c. `R_172_genomes/assemblies.rds` Core database assembly tables.
	 d. `R_172_genomes/seq_region.rds` Core database rank 1, manually filtered
	 seq_region tables used for assessing scaffold N50 values.
   All of these are referenced in downstream code, but only `genome_sizes.txt`
   and `seq_region.rds` were used for analyses included in the publication.
