#!/usr/bin/perl
use warnings;
use Math::Random::MT;
#This script simulates a genotyping by sequencing experiment using genotypes from the inbred or backcross simulator.
#The model of SNP detection now goes from 0 to ntries tries, with p(a1) + p(a2) <= 1.  The rest is p(nothing).
#The most affected loci are expected to be frozen out into superfluous linkage groups or misfit end-caps when Flipper runs on them.
open (DAT, "< $ARGV[0]"); #pagxxivtestafake.txt
open (OUT, "> $ARGV[1]"); #pagxxivtestafakeminhets.txt
open (KLL, "> $ARGV[2]"); #pagxxivtestakilledloci.txt
$datatype = $ARGV[3]; #data_type_f2_self or data_type_backcross
$iseed = $ARGV[4]; #11279929
$nmeth = $ARGV[5]; #400 out of 2000 loci
$methlevel = $ARGV[6]; #0.25
$ntries = $ARGV[7]; #10
if (scalar(@ARGV) != 8) {die "This script needs eight arguments.\n";}
srand($iseed);
for ($i = 0; $i < 624; $i++) { #Set up seeds for the Mersenne Twister.  Each seed is an unsigned integer.
  $v = int(rand(429496));
  for ($j = 0; $j < 4; $j++) {$v .= int(rand(10));}
  $av[$i] = $v;
}
$rng = Math::Random::MT->new(@av);
$currgene = -1;
while ($line = <DAT>) { #Read data into a hash of arrays.
  #The data are assumed to be stacked with a separation of offset between alleles of the same locus.
  chomp $line;
  #>Gene 1216
  #1 1 0 0 1 0 0 0 0 1 1 1 1 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 0 1
  #1 1 0 1 1 1 1 0 1 1 1 0 1 1 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 1 1 0 1 0 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1
  #1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 1 0 1 1 1 0 1 1 1 1 0 0 1 1 0 1 0 1 1 1 1 1 1 0 0 0 0 1 1 1 1 0 1 1 1
  #1 1 1 1 1 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1
  if ($line =~ m/^>/) {
    ($filler, $currgene) = split(/\s+/, $line);
    $deflines{$currgene} = $line;
    $nplants{$currgene} = 0;
  }
  else {
    @vars = split(/ /, $line);
    push(@{$data{$currgene}}, @vars);
    $nplants{$currgene} += scalar(@vars);
  }
}
$datatype =~ s/_/ /g;
$offset = ($currgene + 1) / 2; #The file structure should duplicate every locus for two alleles per locus.
for ($i = 0; $i < $offset; $i++) {$methstatus[$i] = 0;}
$k = 0;
while ($k < $nmeth) {
  $j = int($offset * $rng->rand());
  if ($methstatus[$j] == 0) {
    $methstatus[$j] = 1;
    print KLL "Methylated locus $k = $j\n";
    $k++;
  }
}
close(KLL);
for $key (keys(%nplants)) {
  if (!defined($value)) {$value = $nplants{$key};}
  elsif ($nplants{$key} != $value) {print OUT "Variant number of plants $nplants{$key} was seen for key $key.\n";}
}
print OUT "$datatype\n$nplants{0} $offset 0 ABH\n\n";
#data type f2 self
#200 2000 0 ABH
print "There are $offset loci and $value plants.\n";
for ($i = 0; $i < $value; $i++) {$pconc[$i] = ($rng->rand() + $rng->rand() + $rng->rand()) / 3;}
#for ($i = 0; $i < $value; $i++) {print OUT "pconc[$i] = $pconc[$i]\n";}
$coverage = 0; $nresetsa = 0; $nresetsb = 0;
for $key (sort {$a <=> $b} (keys(%deflines))) {
  if ($key < $offset) {
    for $i (0 .. $#{$data{$key}}) {
      if ($data{$key}[$i] == 1) {
        if ($data{$key+$offset}[$i] == 0) {$pa1 = 1.0; $pa2 = 0;}
        else {$pa1 = 0.5; $pa2 = 0.5;}
      }
      else {
        if ($data{$key+$offset}[$i] == 1) {$pa1 = 0; $pa2 = 1.0;}
        else {
          print "Abending for key = $key data{$key}[$i] = $data{$key}[$i] data{$key+$offset}[$i] = $data{$key+$offset}[$i]\n";
          die "This case of a true null should not occur in these simulations, exiting.\n";
        }
      }
      $pmeth1 = 1;
      if ($methstatus[$key] == 0) {$pmeth2 = 1;}
      else {$pmeth2 = $methlevel;}
      $lim1 = $pconc[$i] * $pa1 * $pmeth1;
      $lim2 = $lim1 + $pconc[$i] * $pa2 * $pmeth2;
      $na = 0; $nb = 0;
      for ($j = 0; $j < $ntries; $j++) {
        $r = $rng->rand();
        if ($r < $lim1) {$na++;}
        elsif ($r < $lim2) {$nb++;}
      }
      $coverage += $na + $nb;
      if ($na == 0) {$data{$key}[$i] = 0; $nresetsa++;}
      if ($nb == 0) {$data{$key+$offset}[$i] = 0; $nresetsb++;} #Otherwise do nothing to this plant for this locus.
    }
  }
}
$nmissing = 0;
for $key (sort {$a <=> $b} (keys(%deflines))) {
  if ($key < $offset) {
    print OUT "*$key\n";
    for $i (0 .. $#{$data{$key}}) {
      if ($data{$key}[$i] == 1) {
        if ($data{$key+$offset}[$i] == 1) {print OUT "H";}
        else {print OUT "A";}
      }
      else {
        if ($data{$key+$offset}[$i] == 1) {print OUT "B";}
        else {print OUT "-"; $nmissing++;} #missing datum when neither allele shows up
        #Note that homozygous B is much more likely to be missed in methylated or differentially amplified loci.
      }
    }
    print OUT "\n";
  }
}
$coverage /= ($offset * $value);
$totalpoints = $offset * $value;
$fracmissing = $nmissing / $totalpoints;
print "The mean read coverage is $coverage x and there are $nmissing missing data values, $fracmissing out of $totalpoints.\n";
