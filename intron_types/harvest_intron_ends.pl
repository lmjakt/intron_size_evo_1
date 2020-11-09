#!/usr/bin/perl -w

## Read files containing intron sequences and extract the first two and the last two introns;
## Then count the number of types of introns that are observed in each species and maybe rank

## The fasta files are in ../family_members/all_genes/
## with one file per species. Each transcript should be represented by a single entry where
## each line of the sequence part represents a single intron sequence. We can thus easily
## get the beginnings and ends of the sequences.

## names of files in: intron_files.txt

$outfile = "Intron_ends_stats.csv";
open($of, ">", $outfile) || die "unable to open $outfile $!";

open(IN, "intron_files.txt") || die "unable to open intron_files.txt for reading $!\n";

while(<IN>){
    chomp;
    push @files, $_;
}




for $f(@files){
    $sp = $f;
    $sp =~ s#.+/([^_]+_[^_]+)_introns.fa#$1#;
    open(IN, $f) || die "unable to open $f for reading $!\n";
    $gene_id = "";
    while(<IN>){
	chomp;
	if($_ =~ /^>(\S+)/){
	    if(length($gene_id)){
		print_entry($of, $sp, $gene_id, $transcript_id, $db, \@lengths, \@ends);
	    }
	    $gene_id = $1;
	    @tmp = split /\t/, $_;
	    $transcript_id = $tmp[1];
	    $db = $tmp[3];
	    @ends = ();
	    @lengths = ();
	    next;
	}
	push @lengths, length($_);
	$e =  substr($_, $lengths[-1]-2, 2);
	push @ends, substr($_, 0, 2).$e;
	$sp_counts{$sp}{$ends[-1]}++;
	$all_counts{$ends[-1]}++;
    }
    print_entry($of, $sp, $gene_id, $transcript_id, $db, \@lengths, \@ends);
}

## print out the counts..
open(OUT, ">", "sp_end_counts.csv") || die "unable to open sp_end_counts.csv $!\n";
for $sp( sort keys %sp_counts ){
    for $end( sort {$sp_counts{$sp}{$a} <=> $sp_counts{$sp}{$b}} keys %{$sp_counts{$sp}} ){
	print OUT "$sp\t$end\t$sp_counts{$sp}{$end}\n";
    }
}

open(OUT, ">", "end_counts.csv") || die "unable to open end_counts.csv $!\n";
for $end( sort {$all_counts{$a} <=> $all_counts{$b}} keys %all_counts ){
    print OUT "$end\t$all_counts{$end}\n";
}

## we can modify the print_entry to also consider rank (from beginning and
## from end.. 
sub print_entry {
    my($of, $sp, $gene_id, $transcript_id, $db, $l_ref, $e_ref) = @_;
    print $of "$sp\t$db\t$gene_id\t$transcript_id\t";
    if(!scalar( @$l_ref)){
	print $of "NA\tNA\n";
    }else{
	print $of join(",", @{$l_ref}), "\t", join(",", @{$e_ref}), "\n";
    }
}
