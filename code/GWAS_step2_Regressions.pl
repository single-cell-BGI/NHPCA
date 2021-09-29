
#!/usr/bin/perl -w
use strict;

open IN, "Traits_list.txt" || die "";
while(<IN>){
    chomp;
    my($a, $b) = split(/\s+/, $_);
    `python ldsc.py --h2-cts $b --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --out celltype/$a --ref-ld-chr-cts celltype.ldcts --w-ld-chr weights_hm3_no_hla/weights.`;
}

