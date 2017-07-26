#!/usr/bin/perl
use warnings;
#This script translates genotypes to a quantitative measure of some trait.
open (DAT, "< $ARGV[0]"); #deduptest_0128b30newfake.txt
open (MEA, "< $ARGV[1]"); #genoquantities.txt
open (OUT, "> $ARGV[2]"); #deduptest_0128b30measurements.txt
while ($line = <MEA>) {
  chomp $line;
  @vars = split(/\s+/, $line);
  push(@actives, $vars[0]);
  push(@falaffects, $vars[1]);
  push(@secaffects, $vars[2]);
}
$currgene = -1;
$parent1measure = 0;
$parent2measure = 0;
for ($i = 0; $i < scalar(@secaffects); $i++) {
  $parent1measure += 2 * $falaffects[$i];
  $parent2measure += 2 * $secaffects[$i];
}
print OUT "Parent 1 measurement = $parent1measure and parent 2 measurement = $parent2measure.\n";
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
$offset = ($currgene + 1) / 2; #The file structure should duplicate every locus for two alleles per locus.
for $key (keys(%nplants)) {
  if (!defined($value)) {$value = $nplants{$key};}
  elsif ($nplants{$key} != $value) {print OUT "Variant number of plants $nplants{$key} was seen for key $key.\n";}
}
print "There are $offset loci and $value plants.\n";
for ($i = 0; $i < $value; $i++) {$measures[$i] = 0;}
for $key (sort {$a <=> $b} (keys(%deflines))) {
  if ($key < $offset) {
    $found = 0;
    for ($j = 0; $j < scalar(@actives); $j++) {
      if ($key eq $actives[$j]) {$found = 1; $savedj = $j;}
    }
    if ($found == 0) {next;}
    for $i (0 .. $#{$data{$key}}) { #i is the plant number.
      if ($data{$key}[$i] == 1) {
        if ($data{$key+$offset}[$i] == 0) {
          if ($found == 1) {$measures[$i] += 2 * $falaffects[$savedj];}
        }
        else {if ($found == 1) {$measures[$i] += $falaffects[$savedj] + $secaffects[$savedj];}}
      }
      else {
        if ($data{$key+$offset}[$i] == 1) {$measures[$i] += 2 * $secaffects[$savedj];}
        else {
          print "Abending for key = $key data{$key}[$i] = $data{$key}[$i] data{$key+$offset}[$i] = $data{$key+$offset}[$i]\n";
          die "This case of a true null should not occur in these simulations, exiting.\n";
        }
      }
    }
  }
}
@sortedmeasures = sort {$a <=> $b} @measures;
for ($i = 0; $i < $value; $i++) {print OUT "$i\t$measures[$i]\t$sortedmeasures[$i]\n";}
