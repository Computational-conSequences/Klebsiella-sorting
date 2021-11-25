#!/usr/bin/env perl

use Getopt::Long;
use strict;

my $speciesTbl = '';
my @groupTbls   = ();
my $resultsDir = 'GroupedSpecies';

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpBrief
    = qq(about:\n)
    . qq(    This program adds a "Same" column with species, genera, or\n)
    . qq(    family labels for genome identifier pairs based on NCBI's\n)
    . qq(    taxonomy\n\n)
    . qq(usage:\n)
    . qq(    $ownName -t <species-table> -g <cluster-files(s)> \n\n)
    . qq(examples:\n)
    . qq(    $ownName -t DataSets/Enterobacteriaceae.species -c Species-Genus/ANI.groups\n\n)
    . qq(options:\n)
    . qq(   -t table with species to genome map [Complete.groups.bz2],\n)
    . qq(       required\n)
    . qq(   -g file(s) with grouped genomes (from pruner.R), required\n)
    . qq(   -o output directory, default '$resultsDir'\n)
    . qq(\n)
    ;

my $options = GetOptions(
    "t=s"    => \$speciesTbl,
    "g=s{,}" => \@groupTbls,
    "o=s"    => \$resultsDir,
) or die $helpBrief;

my $family = $speciesTbl;
$family =~ s{\S+/}{};
$family =~ s{\.\S+}{};

if( !$speciesTbl && !@groupTbls ) {
    die $helpBrief;
}
elsif( !$speciesTbl ) {
    die "I need a genome ID to species table\n$helpBrief";
}
elsif( !$speciesTbl && !@groupTbls ) {
    die "I need pairwise groupance tables\n$helpBrief";
}

my @trueTbls = sort checkFl(uniq(@groupTbls));
my $cFiles = @trueTbls;
if( $cFiles > 0 ) {
    print "working with $cFiles files\n";
}
else {
    die "the files:" . join("\n",uniq(@groupTbls)) ."\ndo not exist\n";
}

print "learning species and genera from $speciesTbl\n";
my( $rSpecies,$rGenus ) = learnTaxa("$speciesTbl");

unless( -d "$resultsDir" ) {
    mkdir("$resultsDir");
}

print "checking groups for unique species\n";
for my $groupTbl ( @trueTbls ) {
    print "   working with $groupTbl\n";
    my $outFile = $groupTbl;
    $outFile =~ s{\S+/}{};
    $outFile =~ s{\.\S+}{};
    $outFile = $resultsDir . "/" . $outFile;
    my( $printed,$rejected )
        = checkGroups($groupTbl,$outFile,$rSpecies,$rGenus);
}

print "        $ownName done\n\n";

#################################################################
#################################################################
####################### sub routines ############################
#################################################################
#################################################################

sub learnTaxa {
    my $spFl = $_[0];
    my %species = ();
    my %genus   = ();
    my $open
        = $spFl =~ m{\.bz2$}    ? "bzip2 -qdc"
        : $spFl =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : "cat";
    open( my $SPFL,"-|","$open $spFl" ) or die "problems with $spFl $!";
  READSPFL:
    while(<$SPFL>) {
        next READSPFL if m{^(Species|#)};
        chomp;
        my( $species,$count,@gcfs ) = split(/\t/,$_);
        my $genus = $species;
        $genus =~ s{^Candidatus\s+}{};
        $genus =~ s{\s+[a-z]+$}{};
        for my $gcf ( @gcfs ) {
            $species{"$gcf"} = $species;
            $genus{"$gcf"}   = $genus;
        }
    }
    close($SPFL);
    return(\%species,\%genus);
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

sub checkFl {
    my @test = @_;
    my $initial = @test;
    my @found   = ();
    my @missing = ();
    for my $groupTbl ( @test ) {
        if( -f $groupTbl ) {
            push(@found,$groupTbl);
        }
        else {
            push(@missing,$groupTbl);
        }
    }
    my $missing = @missing;
    my $found   = @found;
    if( $missing > 0 ) {
        print "$missing missing files:\n",
            join("\n",@missing),"\n";
    }
    if( $found > 0 ) {
        print "found $found files\n";
        return(@found);
    }
    else {
        return();
    }
}

sub naked {
    my $clean = $_[0];
    $clean =~ s{^\S+/}{};
    $clean =~ s{\.\S+$}{};
    return($clean);
}

sub checkGroups {
    my($groupTbl,$outFile,$rSpecies,$rGenus) = @_;
    my $openFl
        = $groupTbl =~ m{\.gz$} ? "gzip -qdc $groupTbl"
        : $groupTbl =~ m{\.bz2$} ? "bzip2 -qdc $groupTbl"
        : "cat $groupTbl";
    my %species = ();
    my %broken  = ();
    my %united  = ();
    open( my $OUTFL,">","$outFile.tmp" );
    print {$OUTFL} join("\t","# Group","Composition"),"\n";
    open( my $INTBL,"-|","$openFl" ) or die "problem with $groupTbl $!";
    while(<$INTBL>) {
        my($groupN,$genomes) = split;
        my @genomes = split(/,/,$genomes);
        my %internal = ();
        my @unknown  = ();
        my $countsp  = 0;
        for my $genome ( @genomes ) {
            if( exists $rSpecies->{"$genome"} ) {
                my $species = $rSpecies->{"$genome"};
                $species{"$species"}++;
                $internal{"$species"}++;
                $countsp++;
            }
            else {
                push(@unknown,$genome);
            }
        }
        my $unknown = @unknown;
        my @species = sort { $internal{"$b"} <=> $internal{"$a"}
                                 || $a cmp $b } keys %internal;
        my $outLine = "";
        for my $species ( @species ) {
            $broken{"$species"}++;
            if( exists $united{"$species"} ) {
                if( $internal{"$species"} > $united{"$species"} ) {
                    $united{"$species"} = $internal{"$species"};
                }
            }
            else {
                $united{"$species"} = $internal{"$species"};
            }
            my $frac = sprintf("%.3f",$internal{"$species"}/$countsp);
            if( length($outLine) > 0 ) {
                $outLine .= "\t" . join("\t",$species,$internal{"$species"},$frac);
            }
            else {
                $outLine .=  join("\t",$species,$internal{"$species"},$frac);
            }
        }
        if( $unknown > 0 ) {
            if( length($outLine) > 0 ) {
                $outLine .= "\t" . join("\t","unknown",$unknown,@unknown);
            }
            else {
                $outLine .= join("\t","unknown",$unknown,@unknown);
            }
        }
        print {$OUTFL} join("\t",$groupN,$outLine),"\n";
    }
    close(INTBL);
    close($OUTFL);
    rename("$outFile.tmp","$outFile.composition");
    open( my $OUTFL2,">","$outFile.tmp" );
    print {$OUTFL2} join("\t","# Name","Count","Groups","Maximum","Proportion"),"\n";
    my @species = sort { $species{"$b"} <=> $species{"$a"}
                             || $a cmp $b } keys %species;
    for my $species ( @species ) {
        my $maxPC = sprintf("%.3f",$united{"$species"}/$species{"$species"});
        print {$OUTFL2} join("\t",
                            $species,$species{"$species"},
                            $broken{"$species"},
                            $united{"$species"},$maxPC),"\n";
    }
    close($OUTFL2);
    rename("$outFile.tmp","$outFile.breaks");
}
