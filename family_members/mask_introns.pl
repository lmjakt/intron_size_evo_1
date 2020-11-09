#!/usr/bin/perl -w
use DBI;
use List::Util qw( min max );
## use the repeat_feature and repeat_consensus tables to create a mask for introns for a given
## species.

($db) = $ARGV[0];

## WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## 
## note that this script is not species agnostic as it excludes chromosome names with ALT in it
## this could be made into an optional variable, but as it is I am just going to leave.

## better to do this for all introns for all canonical transcripts; and ordering by position.
## note that this mean that the introns can get mixed up.

$exon_query = 
    "select a.stable_id as gene, a.seq_region_start as start, a.seq_region_end as end, b.name as chr, a.seq_region_strand as strand, ".
    "c.stable_id as transcript, d.rank, e.stable_id as exon, e.seq_region_start as ex_start, e.seq_region_end as ex_end ".
    "from gene a ".
    "inner join seq_region b on a.seq_region_id=b.seq_region_id and b.name not regexp 'ALT' ".
    "inner join coord_system f on b.coord_system_id=f.coord_system_id and f.rank < 3 ".
    "inner join transcript c on a.canonical_transcript_id=c.transcript_id ".
    "inner join exon_transcript d on a.canonical_transcript_id=d.transcript_id ".
    "inner join exon e on d.exon_id=e.exon_id ".
    "order by b.name, a.stable_id, e.seq_region_start;";

$rep_query =
    "select b.name as chr, a.seq_region_strand as strand, a.seq_region_start as start, a.seq_region_end as end, ".
    "a.repeat_start as r_start, a.repeat_end as r_end, a.score, ".
    "c.repeat_name, c.repeat_class, c.repeat_type ".
    "from repeat_feature a ".
    "inner join seq_region b on a.seq_region_id=b.seq_region_id and b.name not regexp 'ALT' ".
    "inner join repeat_consensus c on a.repeat_consensus_id=c.repeat_consensus_id ".
    "order by b.name, a.seq_region_start;";

print STDERR "getting exons from Ensembl\n";
%exons = get_table($exon_query, $db);
print STDERR "getting repeats from Ensembl\n";
%repeats = get_table($rep_query, $db);

print STDERR "determining overlaps\n";
## define a fkey veriable
%fkey = %{$exons{keys}};
%rkey = %{$repeats{keys}};
$rep_i = 0;
for($ex_i=0; $ex_i < @{$exons{data}}; $ex_i++){
    $gene = $exons{data}[$ex_i][ $fkey{'gene'} ];
    $transcript = $exons{data}[$ex_i][ $fkey{'transcript'} ];
    $strand = $exons{data}[$ex_i][ $fkey{'strand'} ];
    $ex_rank = $exons{data}[$ex_i][ $fkey{'rank'} ];
    $chr = $exons{data}[$ex_i][ $fkey{'chr'} ];
    $ex_end = $exons{data}[$ex_i][ $fkey{'ex_end'} ];
    ## print STDERR "$gene\t$transcript\t$ex_rank\n";
    ## do we have an intron from this gene following this one?
    if( $ex_i < $#{$exons{data}} && $exons{data}[$ex_i + 1][ $fkey{'gene'} ] eq $gene ){
	## we can define intron begins and ends. The intron rank depends on the
	## the strand
	$int_rank = ($strand == 1) ? $ex_rank : $ex_rank - 1;
	$int_start = $ex_end + 1;
	$int_end = $exons{data}[$ex_i + 1][ $fkey{'ex_start'} ] - 1;
	@int_o_start = ();
	@int_o_end = ();
	@int_o_name = ();
	@int_o_score = ();
	$int_coverage = '0'x(1 + $int_end - $int_start);
	## to ensure that we get all repeats that overlap with the intron
	## first make sure that the chromosome is correct.. This is the only place
	## where we are allowd to cross a chromosome boundary;
	while( $repeats{data}[ $rep_i ][ $rkey{'chr'} ] ne $chr ){
	    $rep_i++;
	}
	## whatever we do here avoid crossing a chromosome
	while( $rep_i > 0 && $repeats{data}[ $rep_i-1 ][ $rkey{'chr'} ] eq $chr
	       && $repeats{data}[ $rep_i ][ $rkey{'end'} ] > $int_start ){
	    $rep_i--;
	}
	## then go up until we have gone far enough..
	while( $rep_i < $#{$repeats{data}} && $repeats{data}[$rep_i + 1][ $rkey{'chr'} ] eq $chr
	       && $repeats{data}[$rep_i][ $rkey{'start'} ] < $int_end ){
	    ## determine if we have an overlap and then define the overlapping region of that
	    $o_start = max( $int_start, $repeats{data}[$rep_i][ $rkey{'start'} ] );
	    $o_end = min( $int_end, $repeats{data}[$rep_i][ $rkey{'end'} ] );
	    ## if the repeat is entirely outside of the range then end < start
	    if($o_end > $o_start){
		push @int_o_start, $o_start;
		push @int_o_end, $o_end;
		push @int_o_name, $repeats{data}[ $rep_i ][ $rkey{'repeat_name'} ];
		push @int_o_score, $repeats{data}[ $rep_i ][ $rkey{'score'} ];
		$l = 1 + $o_end - $o_start;
		substr( $int_coverage, $o_start - $int_start, $l ) = "1"x$l;
	    }
	    $rep_i++;
	}
	print $chr, "\t", $strand, "\t", $gene, "\t", $transcript, "\t", $int_rank,
	    "\t", $int_start, "\t", $int_end, "\t";
	if( @int_o_start ){
	    print join(",", @int_o_name), "\t", join(",", @int_o_score), "\t", join(",", @int_o_start), "\t", join(",", @int_o_end);
	}else{
	    print "NA\tNA\tNA\tNA";
	}
	print "\t$int_coverage\n";
    }
}
    
sub get_table {
    my($query, $db_name) = @_;
    my $dbs = "DBI:mysql:database=%s:host=ensembldb.ensembl.org:port=5306";
    my $db = DBI->connect(sprintf($dbs, $db_name), "anonymous", "") || die $DBI::errstr;
    my $sth = $db->prepare($query);
    print STDERR "Executing query:\n$query\n";
    $sth->execute();
    my $names = $sth->{NAME};
    my $numFields = $sth->{'NUM_OF_FIELDS'};
    my %data = ();
    $data{fields} = [@$names];
    for my $i(0..$#{$data{fields}}){
	$data{keys}{ $names->[$i] } = $i;
    }
    print STDERR "getting data\n";
    my $line_no=0;
    while(my $ref = $sth->fetchrow_arrayref){
	for(my $i=0; $i < $numFields; $i++){
	    $data{data}[$line_no][$i] = defined($$ref[$i]) ? $$ref[$i] : "NULL";
	}
	$line_no++;
    }
    $db->disconnect();
    return(%data);
}
