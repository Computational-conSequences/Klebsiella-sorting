#!/usr/bin/env perl

use strict;
use Getopt::Long;

my @groupfiles = ();
my $results    = 'groupComparison.pairs';

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpBrief
    = qq(about:\n)
    . qq(    This program compares the groups resulting from cutting\n)
    . qq(    a hierarchy using the optimal thresholds for different\n)
    . qq(    genome distance measures\n\n)
    . qq(usage:\n)
    . qq(    $ownName -i <group files> -o <results file> \n\n)
    . qq(examples:\n)
    . qq(    $ownName -t Species-Genus/*.groups -o compTables.pairs\n\n)
    . qq(options:\n)
    . qq(   -i group files from prune.R [Species-Genus/*.groups], required\n)
    . qq(   -o output file, will be bzipped, default '$results'\n)
    . qq(\n)
    ;

my $options = GetOptions(
    "i=s{,}" => \@groupfiles,
    "o=s"    => \$results,
) or die $helpBrief;

if( !@groupfiles ) {
    die $helpBrief;
}

print "examining groups\n";
my %pairs = ();
my @names = ();
for my $infile ( @groupfiles ) {
    print "   extracting data from $infile\n";
    my($name,$pairs) = readPairs($infile);
    push(@names,$name);
    $pairs{"$name"} = $pairs;
}

print "figuring out all pairs\n";
my %countall = mixHashes(\@names,\%pairs);

my $outfile = "$results";
print "saving results to $outfile.bz2\n";
open( my $OUTFL,">","$outfile" ) or die "problems with $outfile $!";
print {$OUTFL} join("\t","Pair",@names),"\n";
for my $pair ( sort keys %countall ) {
    my @vector = ();
    for my $name ( @names ) {
        my $add = exists $pairs{"$name"}->{"$pair"} ? 1 : 0;
        push(@vector,$add);
    }
    print {$OUTFL} join("\t",$pair,@vector),"\n";
}
system(qq(bzip2 -f --best $outfile));

print "done with $ownName\n\n";

###########################################################
#################### subroutines ##########################
###########################################################

sub readPairs {
    my $infile = $_[0];
    my $name = $infile;
    $name =~ s{\S+/}{};
    $name =~ s{\.\S+}{};
    my %genericHash = ();
    open( my $INFL,"<","$infile" );
    while(<$INFL>) {
        my($groupid,$group) = split;
        my @group = sort split(/,/,$group);
        while( my $gcf1 = shift @group ) {
            for my $gcf2 ( @group ) {
                my $pair = join(",",$gcf1,$gcf2);
                $genericHash{"$pair"} = 1;
            }
        }
    }
    close($INFL);
    return($name,\%genericHash);
}

sub mixHashes {
    my ($rnames,$rpairs) = @_;
    my %mixed = ();
    for my $name ( @{ $rnames } ) {
        print "   adding pairs from $name\n";
        for my $key ( keys %{ $rpairs->{"$name"} } ) {
            $mixed{"$key"}++;
        }
    }
    return(%mixed);
}
