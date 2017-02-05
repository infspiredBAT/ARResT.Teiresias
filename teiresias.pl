#!/usr/bin/perl -w

use Getopt::Long;

$| = 1;

$version = localtime((stat($0))[9]);

$avatar = "# ARResT/Teiresias - discovering new subsets of stereotyped antigen receptor sequences
# Agathangelidis and Darzentas
# $version
# http://bat.infspire.org | bat\@infspire.org";

$usage = qq($avatar

teiresias.pl

--in		-> table: tab-delimited, NO header line, Unix ASCII - see demo.table
                 1st column:  required!  label, make sure they're unique
                 2nd column:  required!  amino acid junction sequence, with anchors
                 3rd column:  optional   IMGT-formatted V gene (and allele)
                 4th column:  optional   annotation

    other (optional) arguments - defaults are current standards for CLL (PMID:22415752)
--ide       -> minimum identity [0-1], default: 0.5
--sim       -> minimum similarity [0-1], default: 0.7
--lendiff   -> max length difference [0,1,2,3], default: 0
--offdiff   -> max offset difference [0,1,2,3], default: 0
--nophylo   -> ignore phylogenetic relationships as defined in 'phyloequiv'

--demo      -> run demo (~30sec in a multicore system)
--debug     -> do not erase intermediate files


    note / hints / tips
TEIRESIAS and NR-grep are OK on Ubuntu 64bit - otherwise download/compile from original sources
 - be aware though that the latest TEIRESIAS version gives unexpected results.
When possible, sequences in results table are uppercased with the pattern that clustered them.
You can find the same sequence in different clusters (i.e. ~fuzzy clustering)
 - these are used as linkers for higher-level clustering (L1, L2, etc),
   where linkage can actually be very loose so be careful with interpretation.
);

use Cwd 'abs_path';
$path = (split "/teiresias.pl",abs_path($0))[0];
$methoddir = $path;

# arg 0 : aa to ignore
# arg 1 : IBM's TEIRESIAS' l
# arg 2 : IBM's TEIRESIAS' w
# arg 3 : IBM's TEIRESIAS' k
# arg 4 : IBM's TEIRESIAS' n
# arg 5 : IMGT's aa equivalences
# arg 6 : minimum identity   [0-1]
# arg 7 : minimum similarity [0-1]
# arg 8 : max length difference [0,1,2,3]
# arg 9 : max offset difference [0,1,2,3]
# arg 10: score threshold
# arg 11: minimum cluster overlap [0-1], i.e. minimum overlap between clusters to merge
# arg 12: phylo equivalences
# defaults:
$args[0]  = 2;
$args[1]  = 3;
$args[2]  = 6;
$args[3]  = 2;
$args[4]  = 1;
$args[5]  = "$methoddir/imgt11.aaequiv";
$args[6]  = 0.5; # 0.5
$args[7]  = 0.7; # 0.7
$args[8]  = 0; # 0
$args[9]  = 0; # 0
$args[10] = 0;
$args[11] = 0;
$args[12] = "$methoddir/phyloequiv";

my $options_res = GetOptions(
    "in:s"          =>  \$in,
    "ide:s"         =>  \$ide,
    "sim:s"         =>  \$sim,
    "lendiff:s"     =>  \$lendiff,
    "offdiff:s"     =>  \$offdiff,
    "nophylo!"      =>  \$nophylo,
    "demo!"         =>  \$demo,
    "debug!"        =>  \$debug
    );

if ($demo) { $in = "$path/demo.table"; }
     unless ($in) { die $usage; }

if (defined($ide)) {            $args[6]    = $ide; }
if (defined($sim)) {            $args[7]    = $sim; }
if (defined($lendiff)) {        $args[8]    = $lendiff; }
if (defined($offdiff)) {        $args[9]    = $offdiff; }
if ($nophylo) {                 $args[12]   = ''; }
unless (defined($annotation)) { $annotation = ''; }
                                $out_table  = "$in.teiresias.table";
                                # .ide0.5-sim0.7-len0off0-nophylo
                                $out_table_args  = "$out_table.ide$args[6]-sim$args[7]-len$args[8]"."off$args[9]";
                                $out_table_args .= "-nophylo" if $args[12] eq '';

print "$avatar\n\n";


print "(i) dealing with the input\n";
    #open (IN, $in);
    #while (<IN>) {
    #    if (!/^label/) {
    #        chomp;
    #        @data = split(/\t/,$_,-1);
    #        print "$data[0]\t$data[2]\t$data[1]\t\n";
    #    }
    #}
    #close IN;
    #exit;
    chomp($counts{in}      = `wc -l $in | cut -f1 -d' '`);
    chomp($counts{in_uniq} = `cut -f1 $in | sort -u | wc -l`);
    $out = `perl $methoddir/teiresias--input.pl $in $args[0]`;
    if (!-e "$in.ign$args[0].dat" || -s "$in.ign$args[0].dat" == 0) { print "(!) no usable data in table - exiting\n"; exit; }


print "(i) running TEIRESIAS and going through patterns\n";

    if (-e "$in.ign$args[0].thr.hits") { print "   (!) .thr.hits exists, delete if you want to rerun - skipping to clustering\n"; goto SKIP1; }

    $c = $args[1] - 1;
    `$methoddir/teiresias -i$in.ign$args[0].dat -o$in.ign$args[0].thr -l$args[1] -c$c -k$args[3] -w$args[2] -n$args[4] -v -b$args[5]`;
    if (!-e "$in.ign$args[0].thr") { print "(!) no IBM TEIRESIAS output - exiting\n"; exit; } else { chomp($patterns = `grep -c -v '^#' $in.ign$args[0].thr`); if ($patterns == 0) { print "(!) no patterns from IBM's TEIRESIAS, try another dataset - exiting\n"; exit; }}

    $out = `perl $methoddir/teiresias--nrgrep.pl $in.ign$args[0].thr $in.ign$args[0].board4nrgrep $args[6] $args[7] $methoddir $patterns`;
    chomp($flt_patterns = `grep -c '^#' $in.ign$args[0].thr.hits`); if ($flt_patterns == 0) { print STDERR "(!) no patterns left after filtering, maybe try more relaxed parameters - exiting\n"; exit; }

SKIP1:


print "(i) clustering\n";

    $out = `perl $methoddir/teiresias--pairs.pl      $in.ign$args[0].thr.hits $in.orig.board $args[6] $args[7] $args[0] $args[8] $args[9] $args[12]`;
    if (!-e "$in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout" || -s "$in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout" == 0) { print STDERR "(!) no pairs - exiting\n"; exit; }

    $out = `perl $methoddir/teiresias--clustering.pl $in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout $in.orig.board $args[10] $args[11]`;
    if (!-e "$in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout.$args[10].ratio$args[11].lite.clusters.table" || -s "$in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout.$args[10].ratio$args[11].lite.clusters.table" == 0) { print STDERR "(!) no clusters - exiting\n"; exit; }
    chomp($counts{clustered} = `grep -v '^#' $in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout.$args[10].ratio$args[11].lite.clusters.table | cut -f1 | sort -u | wc -l`);
    chomp($counts{clusters}  = `grep -v '^#' $in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout.$args[10].ratio$args[11].lite.clusters.table | cut -f7 | sort -u | wc -l`);


print "(=) done, see $out_table [also copied to $out_table_args]\n";

    $args4table = qq($avatar
# arg 0 : aa to ignore                  $args[0]
# arg 1 : IBM's TEIRESIAS' l            $args[1]
# arg 2 : IBM's TEIRESIAS' w            $args[2]
# arg 3 : IBM's TEIRESIAS' k            $args[3]
# arg 4 : IBM's TEIRESIAS' n            $args[4]
# arg 5 : IMGT's aa equivalences        $args[5]
# arg 6 : identity threshold            $args[6]
# arg 7 : similarity threshold          $args[7]
# arg 8 : length diff                   $args[8]
# arg 9 : offset diff                   $args[9]
# arg 10: score threshold               $args[10]
# arg 11: cluster link ratio            $args[11]
# arg 12: phylo equivalences            $args[12]);

    $counts{clustered100} = sprintf "%.2f", 100 * $counts{clustered} / $counts{in};
    $counts4table = qq(
# in        all     : $counts{in}  ($counts{in_uniq} unique)
# clustered abs | % : $counts{clustered}\t| $counts{clustered100}%
# clusters       L0 : $counts{clusters});
    print "$counts4table\n";

    `echo "$args4table$counts4table" | cat - $in.ign$args[0].thr.hits.ide$args[6]sim$args[7].pairs.layout.$args[10].ratio$args[11].lite.clusters.table > $out_table`;
    `cp -p $out_table $out_table_args`;
    `rm $in.ign* $in.orig.board` unless $debug;


exit;
