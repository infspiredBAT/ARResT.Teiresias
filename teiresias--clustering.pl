#!/usr/bin/perl -w

$| = 1;

use Sort::Naturally;

$cluster_link_ratio = $ARGV[3];

open (BOARD, "$ARGV[1]") || die "Cannot read $ARGV[1]";
while (<BOARD>) {
    if (/^>/) {
        chomp;
        $sample = (split /\t/)[0];
        $sample =~ s/>//;
        $sample =~ s/ //g;
        $samples{$sample} = (split /\t/)[1];
        $genes{$sample} = (split /\t/)[2];
    } elsif (/\w/) {
        chomp;
        $lengths{$sample} = length($_);
    }
}
close BOARD;

$threshold = $ARGV[2];

open (LITECLUSTERS, ">$ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters") || die "Cannot create $ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters";
open (FILE, "$ARGV[0]") || die "Cannot read $ARGV[0]";
%numbers = ();
while (<FILE>) {
    chomp;
    $_ =~ s/\"//g;

    $one = (split "\t")[0];
    $two = (split "\t")[1];
    @one_split = (split "  ", $one, -1);
    @two_split = (split "  ", $two, -1);
    $sample_one = $one_split[0];
    $sample_two = $two_split[0];
    $annotation{$sample_one} = pop @one_split;
    $annotation{$sample_two} = pop @two_split;

    $score = (split "\t")[2];
    if ($score >= $threshold) {
        $per_score{$score}{$sample_one}{$sample_two}{one} = $one;
        $per_score{$score}{$sample_one}{$sample_two}{two} = $two;
        $per_score_single{$score}{$sample_one} = $one;
        $per_score_single{$score}{$sample_two} = $two;
        unless (defined($pairs{$sample_one}{$sample_two})) {
            ++$connectivity{$sample_one};
            ++$connectivity{$sample_two};
        }
        $pairs{$sample_one}{$sample_two} = 1;
        unless (defined($numbers{vertices_mem}{$sample_one})) {
            ++$numbers{vertices_count};
            $numbers{vertices_mem}{$sample_one} = 1;
            if (defined($samples{$sample_one})) {
                ++$numbers{all_vertices_count};
            }
        }
        unless (defined($numbers{vertices_mem}{$sample_two})) {
            ++$numbers{vertices_count};
            $numbers{vertices_mem}{$sample_two} = 2;
            if (defined($samples{$sample_two})) {
                ++$numbers{all_vertices_count};
            }
        }
        ++$numbers{edges};
    }
}
close FILE;

foreach $score (sort {$b<=>$a} keys %per_score) {
    foreach $sample_one (sort keys %{$per_score{$score}}) {
        foreach $sample_two (sort keys %{$per_score{$score}{$sample_one}}) {
            $per_score_extra{$score}{$connectivity{$sample_one}}{$lengths{$sample_one}}{$sample_one}{$connectivity{$sample_two}}{$lengths{$sample_two}}{$sample_two} = 1;
        }
    }
}

$level = 0;

AGAIN:

$clusterindex = -1;

foreach $score (sort {$b<=>$a} keys %per_score_extra) {
    foreach $conne_one (sort {$b<=>$a} keys %{$per_score_extra{$score}}) {
        foreach $len_one (sort {$b<=>$a} keys %{$per_score_extra{$score}{$conne_one}}) {
            foreach $sample_one (keys %{$per_score_extra{$score}{$conne_one}{$len_one}}) {
                foreach $conne_two (sort {$b<=>$a} keys %{$per_score_extra{$score}{$conne_one}{$len_one}{$sample_one}}) {
                    foreach $len_two (sort {$b<=>$a} keys %{$per_score_extra{$score}{$conne_one}{$len_one}{$sample_one}{$conne_two}}) {
                        foreach $sample_two (keys %{$per_score_extra{$score}{$conne_one}{$len_one}{$sample_one}{$conne_two}{$len_two}}) {

                            unless ($sample_one eq $sample_two) {

                                if (!defined($member2cluster{$sample_one}) && !defined($member2cluster{$sample_two})) {

                                    ++$clusterindex;
                                    $cluster2member{$clusterindex}{$score}{$sample_one} = $per_score{$score}{$sample_one}{$sample_two}{one};
                                    ++$clustersize{$clusterindex};
                                    $cluster2member{$clusterindex}{$score}{$sample_two} = $per_score{$score}{$sample_one}{$sample_two}{two};
                                    ++$clustersize{$clusterindex};
                                    $member2cluster{$sample_one}{$score}{$clusterindex} = 1;
                                    $member2cluster{$sample_two}{$score}{$clusterindex} = 1;

                                } elsif (defined($member2cluster{$sample_one}) && !defined($member2cluster{$sample_two})) {

                                    if ($level <= 1) {
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_one}}) {
                                            foreach $onecluster (keys %{$member2cluster{$sample_one}{$score1}}) {
                                                $flag = 1;
                                                foreach $score2 (keys %{$cluster2member{$onecluster}}) {
                                                    foreach $oneclustermember (keys %{$cluster2member{$onecluster}{$score2}}) {
                                                        if (!defined($pairs{$oneclustermember}{$sample_two}) && !defined($pairs{$sample_two}{$oneclustermember})) {
                                                            $flag = 0;
                                                        }
                                                    }
                                                }
                                                if ($flag) {
                                                    $cluster2member{$onecluster}{$score}{$sample_two} = $per_score{$score}{$sample_one}{$sample_two}{two};
                                                    ++$clustersize{$onecluster};
                                                    $member2cluster{$sample_two}{$score}{$onecluster} = 1;
                                                    goto SKIP;
                                                }
                                            }
                                        }
                                        ++$clusterindex;
                                        $cluster2member{$clusterindex}{$score}{$sample_one} = $per_score{$score}{$sample_one}{$sample_two}{one};
                                        ++$clustersize{$clusterindex};
                                        $cluster2member{$clusterindex}{$score}{$sample_two} = $per_score{$score}{$sample_one}{$sample_two}{two};
                                        ++$clustersize{$clusterindex};
                                        $member2cluster{$sample_one}{$score}{$clusterindex} = 1;
                                        $member2cluster{$sample_two}{$score}{$clusterindex} = 1;
                                        
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_one}}) {
                                            foreach $onecluster (keys %{$member2cluster{$sample_one}{$score1}}) {
                                                $cluster2cluster{$onecluster}{$clusterindex}{sum} += $score;
                                                ++$cluster2cluster{$onecluster}{$clusterindex}{count};
                                                ++$connectivity{$onecluster};
                                                ++$connectivity{$clusterindex};
                                                if (!defined($cluster2cluster{$onecluster}{$clusterindex}{members}{$sample_one}) && !defined($cluster2cluster{$onecluster}{$clusterindex}{members}{$sample_two})) {
                                                    ++$cluster2cluster{$onecluster}{$clusterindex}{members_count};
                                                    $cluster2cluster{$onecluster}{$clusterindex}{members}{$sample_one} = 1;
                                                    $cluster2cluster{$onecluster}{$clusterindex}{members}{$sample_two} = 1;
                                                }
                                            }
                                        }
                                      SKIP:
                                    } else {
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_one}}) {
                                            foreach $onecluster (keys %{$member2cluster{$sample_one}{$score1}}) {
                                                $cluster2member{$onecluster}{$score}{$sample_two} = $per_score{$score}{$sample_one}{$sample_two}{two};
                                                ++$clustersize{$onecluster};
                                                $member2cluster{$sample_two}{$score}{$onecluster} = 1;
                                            }
                                        }
                                    }
                                    
                                } elsif (!defined($member2cluster{$sample_one}) && defined($member2cluster{$sample_two})) {
                                    
                                    if ($level <= 1) {
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_two}}) {
                                            foreach $twocluster (keys %{$member2cluster{$sample_two}{$score1}}) {
                                                $flag = 1;
                                                foreach $score2 (keys %{$cluster2member{$twocluster}}) {
                                                    foreach $twoclustermember (keys %{$cluster2member{$twocluster}{$score2}}) {
                                                        if (!defined($pairs{$twoclustermember}{$sample_one}) && !defined($pairs{$sample_one}{$twoclustermember})) {
                                                            $flag = 0;
                                                        }
                                                    }
                                                }
                                                if ($flag) {
                                                    $cluster2member{$twocluster}{$score}{$sample_one} = $per_score{$score}{$sample_one}{$sample_two}{one};
                                                    ++$clustersize{$twocluster};
                                                    $member2cluster{$sample_one}{$score}{$twocluster} = 1;
                                                    goto SKIP;
                                                }
                                            }
                                        }
                                        ++$clusterindex;
                                        $cluster2member{$clusterindex}{$score}{$sample_one} = $per_score{$score}{$sample_one}{$sample_two}{one};
                                        ++$clustersize{$clusterindex};
                                        $cluster2member{$clusterindex}{$score}{$sample_two} = $per_score{$score}{$sample_one}{$sample_two}{two};
                                        ++$clustersize{$clusterindex};
                                        $member2cluster{$sample_one}{$score}{$clusterindex} = 1;
                                        $member2cluster{$sample_two}{$score}{$clusterindex} = 1;
                                        
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_two}}) {
                                            foreach $twocluster (keys %{$member2cluster{$sample_two}{$score1}}) {
                                                $cluster2cluster{$twocluster}{$clusterindex}{sum} += $score;
                                                ++$cluster2cluster{$twocluster}{$clusterindex}{count};
                                                ++$connectivity{$twocluster};
                                                ++$connectivity{$clusterindex};
                                                if (!defined($cluster2cluster{$twocluster}{$clusterindex}{members}{$sample_one}) && !defined($cluster2cluster{$twocluster}{$clusterindex}{members}{$sample_two})) {
                                                    ++$cluster2cluster{$twocluster}{$clusterindex}{members_count};
                                                    $cluster2cluster{$twocluster}{$clusterindex}{members}{$sample_one} = 1;
                                                    $cluster2cluster{$twocluster}{$clusterindex}{members}{$sample_two} = 1;
                                                }
                                            }
                                        }
                                      SKIP:
                                    } else {
                                        foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_two}}) {
                                            foreach $twocluster (keys %{$member2cluster{$sample_two}{$score1}}) {
                                                $cluster2member{$twocluster}{$score}{$sample_one} = $per_score{$score}{$sample_one}{$sample_two}{one};
                                                ++$clustersize{$twocluster};
                                                $member2cluster{$sample_one}{$score}{$twocluster} = 1;
                                            }
                                        }
                                    }
                                    
                                } elsif (defined($member2cluster{$sample_one}) && defined($member2cluster{$sample_two})) {
                                    
                                    foreach $score1 (sort {$b<=>$a} keys %{$member2cluster{$sample_one}}) {
                                        foreach $onecluster (keys %{$member2cluster{$sample_one}{$score1}}) {
                                            foreach $score2 (sort {$b<=>$a} keys %{$member2cluster{$sample_two}}) {
                                                foreach $twocluster (keys %{$member2cluster{$sample_two}{$score2}}) {
                                                    $cluster2cluster{$onecluster}{$twocluster}{sum} += $score;
                                                    ++$cluster2cluster{$onecluster}{$twocluster}{count};
                                                    ++$connectivity{$onecluster};
                                                    ++$connectivity{$twocluster};
                                                    if (!defined($cluster2cluster{$onecluster}{$twocluster}{members}{$sample_one}) && !defined($cluster2cluster{$onecluster}{$twocluster}{members}{$sample_two})) {
                                                        ++$cluster2cluster{$onecluster}{$twocluster}{members_count};
                                                        $cluster2cluster{$onecluster}{$twocluster}{members}{$sample_one} = 1;
                                                        $cluster2cluster{$onecluster}{$twocluster}{members}{$sample_two} = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
                
foreach $cluster (sort {$a<=>$b} keys %cluster2member) {
    if ($clustersize{$cluster} > 1) {
        $cluster2print = "$level-" . sprintf "%04d", $cluster;
        foreach $score (sort {$b<=>$a} keys %{$cluster2member{$cluster}}) {
            foreach $member (sort keys %{$cluster2member{$cluster}{$score}}) {
                if ($level > 0) {
                    $previouslevel = $level - 1;
                    $member2print = "$previouslevel-" . sprintf "%04d", $member;
                } else {
                    $member2print = "\"$cluster2member{$cluster}{$score}{$member}\"";
                }
                $score2print = sprintf "%.0f", $score;
                print LITECLUSTERS "$cluster2print\t$member2print\t$score2print\n";
                if ($level == 0) {
                    if (defined($annotation{$member})) {
                        ++$clusterspec{$cluster2print}{$annotation{$member}};
                    }
                }
            }
        }
    }
}

%clusterspec = ();
%per_score = ();
%per_score_extra = ();

foreach $onecluster (keys %cluster2cluster) {
    foreach $twocluster (keys %{$cluster2cluster{$onecluster}}) {
        if ($cluster2cluster{$onecluster}{$twocluster}{members_count}/$clustersize{$onecluster} >= $cluster_link_ratio && $cluster2cluster{$onecluster}{$twocluster}{members_count}/$clustersize{$twocluster} >= $cluster_link_ratio) {
            $score = sprintf "%.0f", $cluster2cluster{$onecluster}{$twocluster}{sum}/$cluster2cluster{$onecluster}{$twocluster}{count};
            if ($score > $threshold) {
                $per_score{$score}{$onecluster}{$twocluster} = 1;
                $per_score_extra{$score}{$connectivity{$onecluster}}{$clustersize{$onecluster}}{$onecluster}{$connectivity{$twocluster}}{$clustersize{$twocluster}}{$twocluster} = 1;
            }
        } elsif ($level >= 1) {
            $score = sprintf "%.0f", $cluster2cluster{$onecluster}{$twocluster}{sum}/$cluster2cluster{$onecluster}{$twocluster}{count};
            if ($score > $threshold) {
                $per_score{$score}{$onecluster}{$twocluster} = 1;
                $per_score_extra{$score}{$connectivity{$onecluster}}{$clustersize{$onecluster}}{$onecluster}{$connectivity{$twocluster}}{$clustersize{$twocluster}}{$twocluster} = 1;
            }
        }
    }
}

if ($clusterindex >= 1) {
    ++$level;
    %cluster2member = ();
    %member2cluster = ();
    %cluster2cluster = ();
    %clustersize = ();
    goto AGAIN;
}

close LITECLUSTERS;

open (BLO, "$ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters") || die "Cannot read $ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters";
while (<BLO>) {
    chomp;
    $_ =~ s/\"//g;
    if (/^0/) {
        $node = (split /\t/)[1];
        $sample = (split /  /, $node)[0];
        $mod_node = (split /  /, $node)[0] . "\t" . (split /  /, $node)[1] . "\t" . (split /  /, $node)[2] . "\t" . "$annotation{$sample}";
        $cluster = (split /\t/)[0];
        $level = (split /-/, $cluster)[0];
        $excelify{$level}{$mod_node}{$cluster} = (split /\t/)[2];
    } elsif (/[1-9]+-/) {
        $cluster1 = (split /\t/)[0];
        $cluster2 = (split /\t/)[1];
        $level = (split /-/, $cluster1)[0];
        $excelify{$level}{$cluster2}{$cluster1} = (split /\t/)[2];
    }
}
close BLO;
$maxlevel = $level;

foreach $a (sort keys %{$excelify{0}}) {
    %alters = ();
    $alters{"$a\t"} = 1;
    for ($level=0;$level<=$maxlevel;++$level) {
        if (defined($excelify{$level}{$a})) {
            foreach $b (sort keys %{$excelify{$level}{$a}}) {
                foreach $alter (sort keys %alters) {
                    $alter .= "$b\t";
                    $alters_new{$alter} = 1;
                }
                $a = $b;
            }
        } else {
            foreach $alter (sort keys %alters) {
                $alters_new{$alter} = 1;
            }
        }
        %alters = ();
        %alters = %alters_new;
        %alters_new = ();
    }
    foreach $alter (keys %alters) {
        chop $alter;
        $sample = (split /\t/, $alter)[0];
        $sample =~ s/\(/\\(/g;
        $sample =~ s/\)/\\)/g;
        $composcore = (split /\t/, $alter)[2];
        $ide    = substr($composcore,0,6);
        $idepc  = sprintf "%.2f", $ide/10000;
        $sim    = substr($composcore,6,6);
        $simpc  = sprintf "%.2f", $sim/10000;
        $alter  =~ s/$composcore/$idepc\t$simpc\t$genes{$sample}/;
        $l0     = (split /\t/, $alter)[6];
        $alter2order{$l0}{$composcore}{$alter} = 1;
    }
}
open (TABLE, ">$ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters.table") || die "Cannot create $ARGV[0].$threshold.ratio$cluster_link_ratio.lite.clusters.table";
print TABLE "# label\taa junction\tidentity\tsimilarity\t5' gene\tannotation\tL0\tL1\tL2\tL3\tL4\tL5\tetc\n";
foreach $l0 (nsort keys %alter2order) {
    foreach $composcore (sort {$b<=>$a} keys %{$alter2order{$l0}}) {
        foreach $alter (sort keys %{$alter2order{$l0}{$composcore}}) {
            print TABLE $alter . "\n";
        }
    }
}
close TABLE;

exit;
