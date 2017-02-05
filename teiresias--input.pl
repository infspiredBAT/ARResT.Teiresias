#!/usr/bin/perl -w

$| = 1;

$aa2ignore = $ARGV[1];
open (FILE0, $ARGV[0]);
open (AASEQ, ">$ARGV[0].ign$aa2ignore.board");
open (ORIGAASEQ, ">$ARGV[0].orig.board");
while (<FILE0>) {
    #if (/^label/) {
    #    @header = split(/\t/);
    #} elsif (@header && /\w/) {
    if (/^\w/) {
        chomp;
        @data = split(/\t/,$_,-1);
        #for ($i=0;$i<scalar(@header);$i++) { $annotated{$header[$i]} = $data[$i]; }
        #if (   !defined $annotated{'label'}
        #    && !defined $annotated{'V'}
        #    #&& !defined $annotated{'D'}
        #    #&& !defined $annotated{'J'}
        #    #&& !defined $annotated{'V%id'}
        #    && !defined $annotated{'aa junction'}) {
        #    exit;
        #}
        if (   !defined $data[0]
            && !defined $data[1]) {
            print STDERR "(!) first two table columns not defined\n"; exit;
        }
        $annotated{'label'} = $data[0];
        $annotated{'aaju'}  = $data[1];
        if (defined $data[2] && $data[2] =~ /((IG|TR).*)/) {
            $annotated{'5gene'} = $1;
        } else {
            $annotated{'5gene'} = '';
        }
        $annotated{'anno'}  = $data[3] if defined $data[3]; $annotated{'anno'}  = '' if !defined $annotated{'anno'}  || $annotated{'anno'}  !~ /\w/;

        $label = ">$annotated{'label'}\t$annotated{'aaju'}\t$annotated{'5gene'}\t$annotated{'anno'}";

               $aaseq = $annotated{'aaju'};
               $aaseq =~ s/ //g;
        if    ($aaseq =~ /\w/) { print ORIGAASEQ "$label\n$aaseq\n"; }
            $subaaseq = substr $aaseq,$aa2ignore/2,-$aa2ignore/2;
        if ($subaaseq =~ /\w/) { print AASEQ "$label\n$subaaseq\n"; }
    }
}
close ORIGAASEQ;
close AASEQ;

open (FASTA, "$ARGV[0].ign$aa2ignore.board") || die "Cannot read $ARGV[0].ign$aa2ignore.board";
open (DAT,  ">$ARGV[0].ign$aa2ignore.dat") || die "Cannot create $ARGV[0].ign$aa2ignore.dat";
$index = 0;
$check = 0;
while (<FASTA>) {
    chomp;
    if (/^>/) {
        s/^\>+//g;
        s/ +/_/g;
        s/[^A-Za-z0-9]/_/g;
        if ( $check != 0 ) { print DAT "\n"; }	# avoiding printing \n with the first label line			
        print DAT ">$_" . " $index\n";
        ++$index;
        ++$check;
    } else {
        print DAT;
    }
}
print DAT "\n";
close DAT;
close FASTA;

open (FASTA, "$ARGV[0].ign$aa2ignore.board") || die "Cannot read $ARGV[0].ign$aa2ignore.board";
open (BOARD4NRGREP, ">$ARGV[0].ign$aa2ignore.board4nrgrep") || die "Cannot create $ARGV[0].ign$aa2ignore.board4nrgrep";
$check = 0;
while (<FASTA>) {
    chomp;
    if (/^>/) {
        s/[\r\f]//g;
        s/\s.*//;
        if ( $check != 0 ) { print BOARD4NRGREP "\n"; }	# avoiding printing \n with the first label line
        print BOARD4NRGREP "$_" . " \@";
        ++$check;
    } else {
        $_ = uc($_);
        s/[\r\f]//g;			
        s/[^A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X]//g;
        print BOARD4NRGREP $_;
    }
}
print BOARD4NRGREP "\n";
close BOARD4NRGREP;
close FASTA;

exit;
