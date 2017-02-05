#!/usr/bin/perl -w

no warnings 'threads';
use threads;

$| = 1;

$thr_file = $ARGV[0];
$board    = $ARGV[1];
$ide      = $ARGV[2];
$sim      = $ARGV[3];
#$patterns = $ARGV[5];

unless ($thr_file && $board) { die 'oops'; }

$cores   = `nproc`;
#$threads = thread_count($patterns);
$threads = $cores > 2 ? $cores - 2 : $cores;

$i = 0;
open (THR_FILE, "$thr_file") || die "Cannot read $thr_file";
while (<THR_FILE>) {
    unless (/^\#/) {
        $grep_seqlet = (split " ")[2]; # non-eval
        $pattern4len = $grep_seqlet;
        $pattern4len =~ s/(\[[A-Z\.]*\])/x/g;
        $patternlen  = length($pattern4len);
        $identical   = ($pattern4len =~ tr/[A-Z]/i/);
        $equivalent  = ($pattern4len =~ tr/x/e/);
        $similar     = $identical + $equivalent;
        if ($identical/$patternlen >= $ide && $similar/$patternlen >= $sim) {
            $p = $i++ % $threads;
            push @{ $parts[$p] }, "$identical $similar $grep_seqlet";
        }
    }    
}
close THR_FILE;

exit unless $i;

open (BOARD, "$board") || die "Cannot read $board";
while (<BOARD>) {
    chomp;
    if (/^>/) {
        $id = (split " ")[0];
        $seq = (split "\@")[1];
        $lengths{$id} = length($seq);
    }
}
close BOARD;

for ($p=0;$p<$threads;$p++) { push(@$Threads, threads->create(\&nrgrep)); } $hits .= $_->join for @$Threads; undef @$Threads;

open (HITS, ">$thr_file.hits") || die "Cannot create $thr_file.hits";
print HITS $hits;
close HITS;
    
sub nrgrep {
    $phits = '';
    foreach $_ (@{$parts[$p]}) {
        $identical   = (split " ")[0];
        $similar     = (split " ")[1];
        $grep_seqlet = (split " ")[2];
        $grep = `$ARGV[4]/nrgrep -k 1i '@.*$grep_seqlet' $board`;
        $fltgrep = '';
        $fltcount = 0;
        while ($grep =~ /(>\S+)/g) {
            if ($identical/$lengths{$1} >= $ide && $similar/$lengths{$1} >= $sim) {
                $fltgrep .= "$1\n";
                ++$fltcount;
            }
        }
        if ($fltcount > 1) { $phits .= "# $grep_seqlet\n$fltgrep\n"; }
    }
    return $phits;
}

sub thread_count {
    my $e = shift;
       $e >= 1000000 ? 12
     : $e >= 100000  ? 10
     : $e >= 10000   ? 8
     : $e >= 1000    ? 6
     :                 4;
}

exit;
