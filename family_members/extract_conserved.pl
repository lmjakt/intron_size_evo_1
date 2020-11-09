#!/usr/bin/perl -w

use DBI;
use List::Util qw( min max );

## The database is likely to be ensembl_compara_98, but could be a different one
## method_link_species_set_id defines a method used on a set of species
## Suitable identifiers in Ensembl_compara_98 are:
##
# +----------------------------+------------------------------------------------+
# | method_link_species_set_id | name                                           |
# +----------------------------+------------------------------------------------+
# |                       1619 | 27 primates EPO-Low-Coverage                   |
# |                       1621 | 90 eutherian mammals EPO-Low-Coverage          |
# |                       1622 | 90 eutherian mammals GERP Constrained Elements |
# |                       1623 | 33 fish EPO                                    |
# |                       1624 | 60 fish EPO-Low-Coverage                       |
# |                       1625 | 60 fish GERP Constrained Elements              |
# |                       1628 | 13 primates EPO                                |
# |                       1632 | 37 mammals EPO                                 |
# |                       1633 | 53 amniota vertebrates Mercator-Pecan          |
# +----------------------------+------------------------------------------------+
##
## Since these do not have really good names, it is difficult to think of a resonable
## query to identify these.

($db, $method_link_species_set_id) = $ARGV[0];

## to test handling of the data:
$query = "select a.dnafrag_id, a.dnafrag_start, a.dnafrag_end, a.method_link_species_set_id, d.name, d.display_name, c.name, c.length, c.is_reference, b.window_size, b.position, b.diff_score ". 
    "from genomic_align a ".
    "inner join conservation_score b on a.genomic_align_block_id=b.genomic_align_block_id ".
    "inner join dnafrag c on a.dnafrag_id=c.dnafrag_id ".
    "inner join genome_db d on c.genome_db_id=d.genome_db_id limit 100;";

%tmp = get_scores($query, $db);
print join("\t", @{$tmp{fields}}), "\n";
for $i(0..$#{$tmp{data}}){
    for($j=0; $j < $#{$tmp{data}[$i]}; $j++){
	print "\t", $tmp{data}[$i][$j];
    }
    $j = $#{$tmp{data}[$i]};
    $l = (1 + $tmp{data}[$i][2] - $tmp{data}[$i][1] )/$tmp{data}[$i][9];
    @v = unpack( "f".(1 + $l), $tmp{data}[$i][$j] );
    print "\t$l\t", scalar(@v), "\n";
    for($k=0; $k < $l && $k < @v; $k++){ printf("%.2f,", $v[$k]); }
    print "\n";
}


sub get_scores {
    my($query, $db_name) = @_;
    my $dbs = "DBI:mysql:database=%s:host=ensembldb.ensembl.org:port=5306";
##    my $dbs = "DBI:mysql:database=%s:mysql_read_default_file=$ENV{HOME}/.my.cnf";
##    my $db = DBI->connect(sprintf($dbs, $db_name)) || die $DBI::errstr;
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
