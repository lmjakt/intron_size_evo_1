#!/usr/bin/perl -w

## process all_genes/exon_intron_stats.csv
## to extract intron exon coordinates and their respective sizes;

$infile=$ARGV[0];

open(IN, $infile) || die "unable to open $infile $!\n";
while(<IN>){
    chomp;
    @tmp = split;
    ## columns are:
    ## species, database, gene, transcript, chr, strand, gene_coords, exon sizes, intron sizes
    ## coordiantes are nn..nn
    ## big -> small
    ## exons and introns are ordered by transcription so need to be reversed.
    @coords = split /\.\./, $tmp[6];
    $strand = $tmp[5];
    @exon_s = split /,/, $tmp[7];
    if(@exon_s > 1){
	@intron_s = split /,/, $tmp[8];
    }else{
	@intron_s = ();
    }
    if( $strand == -1 ){
	@exon_s = reverse(@exon_s);
	@intron_s = reverse(@intron_s);
    }
    @exon_start = ();
    @intron_start = ();
    $start = $coords[0];
    for $i(0..$#intron_s){
	push @exon_start, $start;
	$start += $exon_s[$i];
	push @intron_start, $start;
	$start += $intron_s[$i]
    }
    push @exon_start, $start if @exon_s > 1;
    print join("\t", @tmp[0..6]), "\t",
	join(",", @exon_start), "\t", join(",", @exon_s), "\t",
	join(",", @intron_start), "\t", join(",", @intron_s),
	"\n";
}
