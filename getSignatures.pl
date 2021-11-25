#!/usr/bin/perl

use strict;
use File::Temp qw( tempfile tempdir ); # to make temporary files/directories
use sigtrap qw(handler signalHandler normal-signals);
use Getopt::Long;
use List::Util qw(sum);

#######################################################
###################### Arguments ######################
#######################################################
my $minSig      = 2;
my $maxSig      = 4;
my $defSig      = 3;
my $defDNA      = 1000;
my $defResDir   = "DNASignature";
my $roundup = "%.4f";
##### controlable
my $fnaFile     = '';
my $signatureLn = $defSig;
my $minDNAln    = $defDNA;
my $resultsDir  = '';
my $circular    = "F";

my $options = GetOptions(
    "i=s" => \$fnaFile,
    "l=s" => \$signatureLn,
    "m=s" => \$minDNAln,
    "o=s" => \$resultsDir,
    "c=s" => \$circular,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if( !$fnaFile ) {
    print qq(about:\n);
    print qq(  This program produces files with DNA oligonucleotide\n)
        . qq(   counts and signatures\n\n);
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i genome file [GCF_000005845.fna.gz], required\n";
    print "   -l signature length, default $defSig\n";
    print "   -m minimum DNA length worth calculating signatures,\n"
        . "       default $defDNA (as per Karlin)\n";
    print "   -c DNA is circular [T/F], default F\n";
    print "   -o output directory, default $defResDir-[signature length]\n";
    print "\n";
    exit;
}

my $fnaRoot = $fnaFile;
$fnaRoot =~ s{^\S*/}{};
$fnaRoot =~ s{\.fna\S*}{};
$fnaRoot =~ s{\.(gz|bz2)}{}g;
my $signatureLn
    = $signatureLn >= $minSig && $signatureLn <= $maxSig ? $signatureLn
    : $defSig;
my $circular = $circular =~ m{^(T|F)$}i ? uc($1) : "F";
my $minDNAln = $minDNAln > 50 ? $minDNAln : $defDNA;

print "\tworking with $fnaRoot signature size -> $signatureLn\n";
print "\t $fnaRoot is circular -> $circular\n";

my $resultsDir
    = length($resultsDir) > 2 ? $resultsDir
    : join("-",$defResDir,$signatureLn);
print "\t saving results at $resultsDir\n";
mkdir($resultsDir) unless(-d $resultsDir);
my $tempFolder = tempdir("/tmp/signatures.XXXXXXXXXXXX");

### pre establish the different oligonucleotides
### (not convenient for oligo sizes above 5)
my @bases = qw( A C G T );
my @nucl  = @bases;
my $current_ln = 0;
while( $current_ln < $signatureLn ) {
    my @grower = ();
    for my $current ( @nucl ) {
        for my $base ( @bases ) {
            my $grown = $current . $base;
            push(@grower,$grown);
        }
    }
    $current_ln = length($grower[0]);
    @nucl = @grower;
}
######### counting both sides:
my %rev_comp = reverse_complement(@nucl);
######### prepared above to count both sides

### learn genome DNA

#### working only with bzipped files
my $calculated = 0;
my $outFile    = "$resultsDir/$fnaRoot.signature.bz2";
my $tmpFile    = "$tempFolder/$fnaRoot.signature.bz2";

my $refSeq = slurpSeq("$fnaFile");

my $nuclHead   = join("\t","Oligonucleotide",@nucl,"Total");
my %gnm_count  = (); ### genomic counts
my $gnmGC      = 0; ### GC content
my $gnmAT      = 0; ### AT content
my $gnm_length = 0;
my @seqIDs
    = sort { length($refSeq->{"$b"}) <=> length($refSeq->{"$a"}) }
    keys %{ $refSeq };
open( my $SIG,"|-","bzip2 -9 >$tmpFile" );
SEQREAD:
for my $seqID ( @seqIDs ) {
    my $seq = $refSeq->{"$seqID"};
    if( length($seq) < $minDNAln ) {
        last SEQREAD;
    }
    print STDERR "\tworking with $seqID\n";
    my($cGC,$cAT,$pGC,$clean) = getGC($seq);
    my($rCounts,$dsCounts)    = countOligos($seq);
    my @calculations
        = calculate_signatures($pGC,$dsCounts);
    $calculated++;
    print {$SIG}
        join("\t","Piece",$seqID,"GC",$cGC,"AT",$cAT,"fGC",$pGC
             ,"Length",length($seq)),"\n";
    print {$SIG} join("\n",$nuclHead,@calculations),"\n";
    ### genomic additions:
    for my $oligo ( keys %{ $dsCounts } ) {
        $gnm_count{"$oligo"} += $dsCounts->{"$oligo"};
    }
    $gnm_length += length($seq);
    $gnmGC      += $cGC;
    $gnmAT      += $cAT;
}

### now the genomic calculations
if( $calculated > 0 ) {
    print STDERR "\tworking the whole thing\n";
    my $pgc = sprintf("$roundup",($gnmGC/($gnmGC + $gnmAT)));
    print {$SIG}
        join("\t","Whole",$fnaRoot,"GC",$gnmGC,"AT",$gnmAT,"fGC",$pgc
             ,"Length",$gnm_length),"\n";
    my @calculations
        = calculate_signatures($pgc,\%gnm_count);
    print {$SIG} join("\n",$nuclHead,@calculations),"\n";
    close($SIG);
    system("mv $tmpFile $outFile 2>/dev/null");
}
else {
    close($SIG);
    unlink("$resultsDir/$fnaRoot.signature.bz2");
}

if( -d "$tempFolder" ) {
    print "\tcleaning up ...\n";
    system "rm -r $tempFolder";
    print "\tdone!\n\n";
}
else {
    print "\ttemp files cleared out\n";
    print "\tdone!\n\n";
}

sub reverse_complement {
    my @nucl = @_;
    my %rev_comp = ();
    for my $nucl ( @nucl ) {
        my $rev_comp = reverse($nucl);
        $rev_comp =~ tr{ATGC}{TACG};
        $rev_comp{"$nucl"} = $rev_comp;
    }
    return(%rev_comp);
}

sub calculate_signatures {
    my($pGC,$nucl_count) = @_;
    my $pAT = 1 - $pGC;
    my $item_total = sum( values %{ $nucl_count } );
    my $total_prop = 0;
    my $total_prob = 0;
    my $total_expt = 0;
    my @line_prop  = ("Proportion");
    my @line_prob  = ("Probability");
    my @line_count = ("Count");
    my @line_expt  = ("Expected");
    my @line_sign  = ("Signature");
    for my $nucl ( @nucl ) {
        my $count = $nucl_count->{"$nucl"};
        push(@line_count,$count);
        ##### now proportions and signatures
        my $prop = $count/$item_total;
        $total_prop += $prop;
        my $print_prop = sprintf("%.5f",$prop);
        push(@line_prop,$print_prop);
        my $gc = $nucl =~ tr{GC}{GC};
        my $at = $nucl =~ tr{AT}{AT};
        #print "$nucl\t( $pGC**$gc ) * ( $pAT**$at )\n";
        my $probability = ( ($pGC/2)**$gc ) * ( ($pAT/2)**$at );
        $total_prob += $probability;
        my $print_prob = sprintf("%.5f",$probability);
        push(@line_prob,$print_prob);
        my $expt_calc = $item_total * $probability;
        my $expt = sprintf("%.1f",$expt_calc);
        push(@line_expt,$expt);
        $total_expt += $expt;
        my $sign = $expt > 0 ? sprintf("%.5f",($count/$expt)) : 0;
        #my $sign = sprintf("%.5f",($prop/$probability));
        push(@line_sign,$sign);
    }
    my @lines_to_print = ();
    my $line_prop  = join("\t",@line_prop,$total_prop);
    my $line_prob  = join("\t",@line_prob,$total_prob);
    my $line_count = join("\t",@line_count,$item_total);
    my $line_expt  = join("\t",@line_expt,$total_expt);
    my $line_sign  = join("\t",@line_sign);
    my @lines_to_print
        = ($line_prop,$line_prob,$line_count,$line_expt,$line_sign);
    return(@lines_to_print);
}

sub signalHandler {
    if( -d "$tempFolder" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $tempFolder";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n";
        die  "    done!\n\n";
    }
}

sub figureCompression {
    my $rootName = $_[0];
    $rootName =~ s{\.(gz|bz2|Z)$}{};
    my $fullName
        = ( -f "$rootName.gz" )  ? "$rootName.gz"
        : ( -f "$rootName.Z" )   ? "$rootName.Z"
        : ( -f "$rootName.bz2" ) ? "$rootName.bz2"
        : ( -f "$rootName" )     ? "$rootName"
        : "none";
    my $catProg
        = $fullName =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : $fullName =~ m{\.bz2$}    ? "bzip2 -qdc"
        : "cat";
    if( $fullName eq "none" ) {
        return();
    }
    else {
        return("$catProg","$fullName");
    }
}

sub slurpSeq {
    my $fnaFl = $_[0];
    my($cat,$truefile) = figureCompression("$fnaFl");
    my %seq = ();
    my $id  = '';
    open( my $FNA,"-|","$cat $truefile" );
    while(<$FNA>) {
        if( m{^>\s*(\S+)} ) {
            $id = $1;
        }
        else {
            chomp;
            $seq{"$id"} .= uc($_);
        }
    }
    close($FNA);
    my $cSeqs = keys %seq;
    if( $cSeqs > 0 ) {
        return(\%seq);
    }
    else {
        return();
    }
}

sub getGC {
    my $seq = $_[0];
    my $countGC = $seq =~ tr{GC}{GC}; ### GC content
    my $countAT = $seq =~ tr{AT}{AT}; ### AT content
    my $fGC     = sprintf("$roundup",($countGC/($countGC + $countAT)));
    my $clean   = length($seq) - ($countGC + $countAT);
    return($countGC,$countAT,$fGC,$clean);
}

sub countOligos {
    my $seq = $_[0];
    my %count_oligo = ();
    my @bases = split(//,$seq);
    ##### procedure to work with circular DNA pieces
    if( $circular eq "T" ) {
        my $add2circle = $signatureLn - 2;
        my @add2circle = @bases[0 .. $add2circle];
        push(@bases,@add2circle);
        print "\tappended: ",join(", ",@add2circle),"\n";
    }
    my $addRange   = $signatureLn - 1;
  BASE:
    for my $n ( 0 .. $#bases ) {
        my $complement = $n + $addRange;
        last if( $complement > $#bases  );
        my $nucl = join "",@bases[$n .. $complement];
        next BASE unless( $nucl =~ m{^[ATGC]+$} );
        #print $nucl,"<-nucleotide\n";exit if($n > 5);
        $count_oligo{"$nucl"}++;
    }
    ####### double strand:
    my %ds_count = ();
    for my $nucl ( @nucl ) {
        my $count_o
            = $count_oligo{"$nucl"} > 0 ? $count_oligo{"$nucl"} : 0;
        my $rev_comp = $rev_comp{"$nucl"};
        my $count = $rev_comp eq $nucl ? $count_o * 2
            : $count_o + $count_oligo{"$rev_comp"};
        $ds_count{"$nucl"}   = $count;
    }
    my $counts = keys %count_oligo;
    if( $counts > 0 ) {
        return(\%count_oligo,\%ds_count);
    }
    else {
        return();
    }
}
