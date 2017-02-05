#!/usr/bin/perl -w

$| = 1;

$ide       = $ARGV[2];
$sim       = $ARGV[3];
$aa2ignore = $ARGV[4];

$sample_count = 0;
open (FILE2, $ARGV[1]);
while (<FILE2>) {
    chomp;
    if (/^>/) {
        $label = $_;
        $label =~ s/^>//;
        @labelarray = split("\t", $label, -1);
        $sample = $labelarray[0];
        $sample =~ s/ +//g;
        $mem{$sample}{label}     = $label;
        $mem{$sample}{vsubgroup} = (split /-/,  $labelarray[2])[0];
        $mem{$sample}{vgene}     = (split /\*/, $labelarray[2])[0];
        ++$sample_count;
    } else {
        $mem{$sample}{aaseq} = lc($_);
        $mem{$sample}{len}   = length($mem{$sample}{aaseq});
    }
}
close FILE2;

if (defined($ARGV[7])) {
    $phylogroup = 1;
    open (PHYLO, $ARGV[7]);
    while (<PHYLO>) {
        chomp;
        @bits = split(/\s+/);
        foreach $bit (@bits) {
            $phyloequivs{$bit} = $phylogroup;
            if ($bit =~ /-/) {       $phylolevel = 'vgene';
            } elsif ($bit =~ /\*/) { $phylolevel = 'vallele';
            } else {                 $phylolevel = 'vsubgroup'; }
        }
        ++$phylogroup;
    }
    close PHYLO;
    $phyloflag = 1;
} else {
    $phyloflag = 0;
}

open (FILE0, $ARGV[0]);
$flag = 0;
$maxpatternscore = 0;
$minpatternscore = 999999999999;
$patternindex = 0;
while (<FILE0>) {
    chomp;
    if (/^\#/ || eof(FILE0)) {
        $patternflag = 0;
        if ($flag && $count > 1) {
            $group = 0;
            foreach $length (sort {$a<=>$b} %per_pattern) {
                foreach $sample (sort keys %{$per_pattern{$length}}) {
                    $groups{$group}{$sample} = 1;
                    ++$finalcount{$group};
                    $grouplensum{$group} += $mem{$sample}{len} - $aa2ignore;
                    $toadd2print{$per_pattern{$mem{$sample}{len}}{$sample}} = 1;
                    $samples_for_blo{$group}{$sample} = 1;
                    $basesample = $sample;
                    for ($i=$mem{$basesample}{len};$i<=$mem{$basesample}{len}+$ARGV[5];++$i) {
                        foreach $sample (sort keys %{$per_pattern{$i}}) {
                            unless ($sample eq $basesample) {
                                if (abs($offsetmem{$basesample}-$offsetmem{$sample}) <= $ARGV[6]) {
                                    if ($similar/($mem{$basesample}{len}-$aa2ignore) >= $sim && $similar/($mem{$sample}{len}-$aa2ignore) >= $sim && $identical/($mem{$basesample}{len}-$aa2ignore) >= $ide && $identical/($mem{$sample}{len}-$aa2ignore) >= $ide) {
                                        $patternblast{$basesample}{$sample} = 1;
                                        $groups{$group}{$sample} = 1;
                                        ++$finalcount{$group};
                                        $grouplensum{$group} += $mem{$sample}{len} - $aa2ignore;
                                        $toadd2print{$per_pattern{$mem{$sample}{len}}{$sample}} = 1;
                                        $samples_for_blo{$group}{$sample} = 1;
                                    }
                                }
                            }
                        }
                    }
                    foreach $line (sort keys %toadd2print) {
                        $toprint{$group} .= $line;
                    }
                    %toadd2print = ();
                    ++$group;
                }
            }
            foreach $group (sort {$a<=>$b} keys %samples_for_blo) {
                if ($finalcount{$group} > 1) {
                    unless (defined($printed{$toprint{$group}})) {
                        $printed{$toprint{$group}} = 1;
                        $toprint{$group} =~ s/\@/$group\t$finalcount{$group}/g;
                        $lenave = $grouplensum{$group} / $finalcount{$group};
                        $ide2print = sprintf "%.4f", $identical/$lenave * 100;
                        $sim2print = sprintf "%.4f", $similar  /$lenave * 100;
#####################################################################
                        if ($ide2print == 100) { $ide4score = '999999';
                        } else {                 $ide4score = $ide2print; }
                        if ($sim2print == 100) { $sim4score = '999999';
                        } else {                 $sim4score = $sim2print; }
                        $patternscore = "$ide4score$sim4score";
                        $patternscore =~ s/\.//g;
#####################################################################
                        $patternscore2print = $patternscore;
                        $toprint{$group} =~ s/\%/$patternscore2print/g;
                        $patternflag = 1;
                        if ($patternscore >= $maxpatternscore) { $maxpatternscore = $patternscore; }
                        if ($patternscore <= $minpatternscore) { $minpatternscore = $patternscore; }
                        foreach $sample (sort keys %{$samples_for_blo{$group}}) {
                            push(@{$samples_only_blo{$patternscore}{$pattern}{$group}},$sample . "  $highlighted{$sample}  $patternscore2print  $annotation{$sample}");
                        }
                    }		    
                }
            }
            if ($patternflag) {
                $patternindices{$pattern} = $patternindex;
                ++$patternindex;
            }
        }
        goto SKIP if eof(FILE0);
        $patternflag = 0;
        $pattern = (split " ")[1];
        $lc_pattern = lc($pattern);
        $pattern4len = $pattern;
        $pattern4len =~ s/(\[[A-Z\.]*\])/x/g;
        $patternlen = length($pattern4len);
        $#patternarray = -1;
        @patternarray = split(//,$pattern4len);
        $identical = ($pattern4len =~ tr/[A-Z]/i/);
        $equivalent = ($pattern4len =~ tr/x/e/);
        $similar = $identical + $equivalent;
        if ($identical/$patternlen >= $ide && $similar/$patternlen >= $sim) {
            $patternflag = 1;
        }
        $flag = 1;
        $count = 0;
        %per_pattern = ();
        %offsetmem = ();
        %finalcount = ();
        %toprint = ();
        %printed = ();
        %samples_for_blo = ();
        %grouplensum = ();
        %groups = ();
        %highlighted = ();
    }

    if (/^>/ && $patternflag) {
        $sample = $_;
        $sample =~ s/^>//;
        $sample =~ s/ +//g;
        $hit = $mem{$sample}{label};
        $aaseq = $mem{$sample}{aaseq};
        @array = split(//,$aaseq);
        while ($aaseq =~ /$lc_pattern/g) {
            $offset = pos($aaseq) - $patternlen + 1;
            $endcheck = length($mem{$sample}{aaseq}) - pos($aaseq) + 1;
            if ($offset > $aa2ignore/2 && $endcheck > $aa2ignore/2) {
                $j = 0;
                for ($i=0;$i<@patternarray;$i++) {
                    if ($patternarray[$i] =~ /[A-Za-z]/) {
                        $array[$offset+$i-1] = uc($array[$offset+$i-1]);
                    }
                }
            } else {
                goto SKIP;
            }
        }
        $aaseq = join("",@array);
        $highlighted{$sample} = $aaseq;
        unless (lc($aaseq) eq $mem{$sample}{aaseq}) { print STDERR "(!) lc(aaseq) ne mem{$sample}{aaseq}\n"; } 
        ++$count;
        $offsetmem{$sample} = $offset;
        $annotation{$sample} = (split "\t", $hit, -1)[3];
        $per_pattern{$mem{$sample}{len}}{$sample} = "$patternindex\t\@\t\%\t$pattern\t$hit\t$aaseq\n";
    }
  SKIP:
}
close FILE0;

%per_pattern = ();

foreach $patternscore (sort {$b <=> $a} keys %samples_only_blo) {
    $normpatternscore = $patternscore;
    foreach $pattern (sort keys %{$samples_only_blo{$patternscore}}) {
        foreach $group (sort {$a<=>$b} keys %{$samples_only_blo{$patternscore}{$pattern}}) {
            for ($i=0;$i<=@{$samples_only_blo{$patternscore}{$pattern}{$group}}-1;++$i) {
                $sample_a = (split "  ", $samples_only_blo{$patternscore}{$pattern}{$group}[$i])[0];
                for ($j=0;$j<=@{$samples_only_blo{$patternscore}{$pattern}{$group}}-1;++$j) {
                    $sample_b = (split "  ", $samples_only_blo{$patternscore}{$pattern}{$group}[$j])[0];
                    
                    $proceed = 1;
                    if ($phyloflag) {
                        if (defined $phyloequivs{$mem{$sample_a}{$phylolevel}} && defined $phyloequivs{$mem{$sample_b}{$phylolevel}}
                                 && $phyloequivs{$mem{$sample_a}{$phylolevel}}     !=     $phyloequivs{$mem{$sample_b}{$phylolevel}}) {
                            $proceed = 0;
                        }
                    }
                    if ($proceed) {
                        if (defined($patternblast{$sample_a}{$sample_b}) || defined($patternblast{$sample_b}{$sample_a})) {
                            unless ($samples_only_blo{$patternscore}{$pattern}{$group}[$i] eq $samples_only_blo{$patternscore}{$pattern}{$group}[$j]) {
                                unless (defined($pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{max}{score})) {
                                    $pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{max}{score} = -1;
                                }
                                if ($normpatternscore >= $pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{max}{score}) {
                                    $pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{max}{score} = $normpatternscore;
                                    $pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{max}{pattern}{$pattern} = 1;
                                }
                                $pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{sum} += $normpatternscore;
                                ++$pairwise_samples{$samples_only_blo{$patternscore}{$pattern}{$group}[$i]}{$samples_only_blo{$patternscore}{$pattern}{$group}[$j]}{count};
                            }
                        }
                    }
                }
            }
        }
    }
}

%patternblast = ();
%samples_only_blo = ();
%mem = ();

$filename = "$ARGV[0].ide$ide" . "sim$sim.pairs.layout";
open (BLO3, ">$filename");
foreach $sample_a (sort keys %pairwise_samples) {
    foreach $sample_b (sort keys %{$pairwise_samples{$sample_a}}) {
        unless (defined($mem{$sample_b}{$sample_a})) {
            $score = $pairwise_samples{$sample_a}{$sample_b}{max}{score};
            $max_patterns = '';
            foreach $pattern (sort keys %{$pairwise_samples{$sample_a}{$sample_b}{max}{pattern}}) {
                $max_patterns .= "$pattern \[$patternindices{$pattern}\],";
            }
            chop $max_patterns;
            print BLO3 "\"$sample_a\"\t\"$sample_b\"\t$score\t$max_patterns\n";
            $mem{$sample_a}{$sample_b} = 1;
        }
    }
}
close BLO3;

exit;
