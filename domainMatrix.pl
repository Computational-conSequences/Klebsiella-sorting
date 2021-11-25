#!/usr/bin/env perl
use strict;

my @fams = qw(
                 pfam
                 cdd
                 cog
                 cd
                 tigrfam
         );

my $famMatch = join("|",@fams);

my $workFam = lc($ARGV[0])
    or die "I need a family to work with [$famMatch]\n\n";
my $workDir = uc("$workFam");
my $outDir  = "famMatrices";
unless( -d "$outDir" ) {
    mkdir("$outDir");
}
opendir( my $FAMD,"$workDir" );
my @files = grep { m{\.$workFam\.} } readdir($FAMD);
closedir($FAMD);

print "learning all families in files\n";
my %totalFam = ();
for my $file ( @files ) {
    my $file2open = "$workDir/$file";
    my %counts = countFams("$file2open");
    for my $fam ( keys %counts ) {
        $totalFam{"$fam"} += $counts{"$fam"};
    }
}
my @families = sort keys %totalFam;

print "building matrix\n";
open( my $MAT,"|-","bzip2 --best > $outDir/$workFam.matrix.bz2" );
print {$MAT} join("\t",@families),"\n";
for my $file ( @files ) {
    my $file2open = "$workDir/$file";
    my %counts = countFams("$file2open");
    my $genome = nakedName("$file");
    my @line = ();
    for my $fam ( @families ) {
        my $current = exists $counts{"$fam"} ? $counts{"$fam"} : 0;
        push(@line,$current);
    }
    print {$MAT} join("\t",$genome,@line),"\n";
}
close($MAT);

sub nakedName {
    my $inName  = $_[0];
    my $outName = $inName;
    $outName =~ s{^\S+/}{};
    $outName =~ s{\.\S+}{};
    if( length($outName) > 0 ) {
        return("$outName");
    }
    else {
        print "there's a problem with $inName -> $outName\n";
        return();
    }
}

sub countFams {
    my $infile = $_[0];
    my %countFam = ();
    open( my $FILE,"-|","bzip2 -qdc $infile" );
    while(<$FILE>) {
        next if( m{^#} );
        my($prots,$fams) = split;
        my @fams = cleanFams("$fams");
        for my $fam ( @fams ) {
            $countFam{"$fam"}++;
        }
    }
    close($FILE);
    return(%countFam);
}

sub cleanFams {
    my $unclean = $_[0];
    my @clean = ();
    for my $toclean ( split(/;/,$unclean) ) {
        $toclean =~ s{\(\S+}{};
        push(@clean,$toclean);
    }
    return(@clean);
}
