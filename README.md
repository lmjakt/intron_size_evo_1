# Scripts and data files for 'Intron size minimisation in teleosts'

These analyses were carried out using a combination of Perl, R and R extensions
written in C. To complicate matters the analyses were carried out on two
separate computers due to operational reasons. An initial R-session was later
split into multiple more clearly defined sessions, but using data objects
exported from the initial session in order to avoid redoing the major parts of
the analysis. Nevertheless, there is an unfortunate amount of overlap between
the initial file and the subsequent targetted analyses.

The files in this git repository come from data extraction and analyses which
were carried out on computer 1 (ws 1). Computer 1 contained the Ensembl databases but
did not have sufficient memory or compute resources for the stages of analysis
depending on large scale sequence alignments.


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


## Source files

- `R_172_genomes/`
  - `R_i72_genomes.R`  
  Initial analyses carried out relying on locally installed databases.
  Provided `genome_sizes.txt`, `seq_region.rds`, `species_lineages.rds` for
  further analyses on ws2.
- `R_common_functions/`
  - `functions.R`  
  Primarily functions for accessing gene coordinate data from Ensembl core
  databases.
- `core_dbs/`  
  Scripts for installing Ensembl databases.
  - `get_databases.sh`  
  Download of Ensembl v. 97 core databases for a preliminary analysis.
  - `uncompress_db_file.sh`  
  Decompression of databases downloaded with `get_database.sh`
  - `install_dbs.sh`  
  Installation of Ensembl v. 97 core databases.
  - `install_db.sh`  
  Installation of a single database after download.
  - `install_ensembl_core_98.sh`  
  Download and installation of Ensembl v. 98 core databases. Calls `install_db.sh`.
- `family_members/`  
  Scripts related to defining a core set of orthologous genes from which to
  define an intron orthology.
  - `compile_all_seqs.pl`  
  - `compile_all_seqs.sh`  
  Extraction of intron and exon sequences from all genes in all
  species. Output two fasta files containing exon and intron sequences respectively.
  - `compile_seqs.pl`  
  - `compile_seqs.sh`  
  Extraction of intron and exon sequences from members of xxxx protein
  families selected in `family_members.R`. Output two fasta files containing
  exons and introns respectively for each family.
  - `count_nucleotides.pl`  
  Determine nucleotide frequencies in intron sequences.
  - `extract_conserved.pl`  
  Identify conserved regions from Ensembl databases. Not used due to unknown
  binary format encoding conservation information.
  - `extract_intron_pos.pl`
  - `extract_intron_pos.sh`  
  Extraction of genomic coordinates, used by
  `R_genome_pos/genome_pos.R`.
  - `family_members.R`  
  Definition of a suitable set of protein families upon which to base the
  intron orthology.
  - `intron_size_exploration.R`  
  Initial analyses for a smaller subset of species. Not included here. Should
  be removed.
  - `mask_introns.pl`  
  - `mask_dr_introns.sh`  
  Create a mask for repetitive regions in *D. rerio*. Analyses not shown as
  the results are somewhat redundant with published reports. Might include in
  later work though.
- `genome_sequences/`  
  Scripts used for initial analyses to obtain version 97 genomes from
  Ensembl. Limited to 15 species.
  - `get_genomes.sh`
  - `get_genomes_2.sh`
- `genome_sequences_98/`  
  Script used to download Ensembl version 98 sequences. Reads in database names
  from output file created in `family_mebers/family_members.R`.
  - `get_genomes_98.sh`
- `intron_types/`  
  Incomplete analyses of sizes of introns stratified by donor and acceptor sites.
  - `harvest_intron_ends.pl`
  - `intron_types.R`
  - `strict_split/`  
    Compiled modified strsplit function.
    - `strict_split.c`
  - `test.R`
    Test compiled `strct_split.c` code.

## Data files

- `R_172_genomes/`  
  Basic data files containing exon and intron coordinates as well as genome
  and assembly properties.
  - `assemblies.rds`  
    Assembly information. Not used for downstream analyses since too variable
    between species.
  - `seq_region.rds`  
    Assembly statistics inferred from `seq_region` tables.
  - `genome_sizes.txt`  
    A list of genome sizes.
  - `species_lineages.rds`  
    A species lineage tree defined from Ensembl compara.
- `family_members/`  
  Intron and exon coordinates and sequences.
  - `all_genes/`
    Intron and exon sequences from all species and all genes. Two sequence
    files for each species containing exons and introns respectively. Exon and
    intron borders are delimited by newlines.
    - `exon_intron_stats.csv`  
	Coordinates of introns and exons.  
	Intron and exon sequences from all species not included in this repository
    but these are available on request. These take the form:
    - `acanthochromis_polyacanthus_exons.fa`
    - `acanthochromis_polyacanthus_introns.fa`
    - `ailuropoda_melanoleuca_exons.fa`
    - `ailuropoda_melanoleuca_introns.fa`
	- etc, with two files per species.
  - `danio_rerio_introns_mask.txt`  
    A repeat mask for introns in *D. rerio* in a tab delimited format. One
    line per intron, with the following columns: chromosome, strand, gene,
    transcript, rank, intron start, intron end, repeat name, repeat score,
    repeat start, repeat end, mask.  
	All repeat columns are comma separated elements with one entry for each
    repeat reported in the intron. The mask is a sequence of 0s ad 1s, where 1
    indicates the presence of a repeat.
  - `db_list.txt`  
  A list of the core Ensembl databases used for the 172 species analyses.
  - `exon_intron_pos.csv`  
  A simplified table of positions and lengths of introns and exons from all
  genes and species. Not used so should be removed from text.
  - `introns_nuc_counts.txt`  
  Nucleotide frequencies. Presumably from intron sequences, but I do not have
  a record of the call. This file has been used to generate random sequences.
  - `orthologue_transcripts/`  
    Exons and intron sequences from all species and 6114 protein
    families. Sequences are organised by family with two files giving exon and
    intron sequences per protein family. Sequence files are not included in
    this archive.
    - `exon_intron_stats.csv`  
	Exon and intron positions and lengths for members of the selected families.
  - `vertebrate_family_members_1_130.txt`  
    Gene identifiers for each species for members of the 6114 protein families
    used for the intron orthology used here.
- `genome_sequences/`  
  Genome sequences obtained for a preliminary analysis. Not included in this repository.
- `intron_types/`  
  Output from an inspection of intron lengths discretized by donor and
  acceptor sequences.
  - `Intron_ends_stats.csv`
  - `end_counts.csv`
  - `intron_files.txt`
  - `sp_end_counts.csv`
