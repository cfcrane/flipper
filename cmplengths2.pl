#!/usr/bin/perl
use warnings;
open(MAP, "< $ARGV[0]"); #deduptest_0219c30fakeminhetsmap.txt
open(ACT, "< $ARGV[1]"); #deduptest_0219c30newouta.txt
open(OUT, "> $ARGV[2]"); #deduptest_0219c30suborder.txt
$searchterm = $ARGV[3]; #forced
while ($line = <ACT>) {
  if ($line !~ m/TELO/) {
    chomp $line;
    @vars = split(/\s+/, $line);
    push(@actuals, $vars[1]);
  }
}
$nloci = scalar(@actuals);
print "Actual locus count = $nloci\n";
$getflag = 0;
while ($line = <MAP>) {
  #5662    5662    -1      -       5.460947        0.010000        0.540000        1.300000
  #linkage group 2 start = 3219    end     1834    count = 248     length = 5.399571
  #3219    3219    4489    0.045455        0.000000        0.045455        0.525253        1.269841
  if ($line =~ m/^\n/) {$getflag = 0;}
  if ($getflag == 1) {
    if ($line =~ m/linkage group/) {
      chomp $line;
      @tars = split(/\s+/, $line);
      $currlg = $tars[2];
      @lg = ();
      for ($i = 0; $i < scalar(@actuals); $i++) {$reals[$i] = -1;}
    }
    else {
      chomp $line;
      @vars = split(/\s+/, $line);
      push(@lg, $vars[1]);
      if ($vars[2] == -1) {
        #Process and output here.
        @reouts = ();
        for ($i = 0; $i < scalar(@actuals); $i++) {
          for ($j = 0; $j < scalar(@lg); $j++) {
            if ($actuals[$i] == $lg[$j]) {push(@reouts, $actuals[$i]); last;}
          }
        }
        if (scalar(@lg) != scalar(@reouts)) {die "lg and reouts arrays differ in size\n";}
        else {
          #The obtained map has 38 linkage groups.
          #For mapped group 0, map start = 4696 and actual start = 2821
          print OUT "For mapped group $currlg\n";
          for ($i = 0; $i < scalar(@lg); $i++) {
            print OUT "$lg[$i]\t$reouts[$i]\n";
          }
        }
      }
    }
  }
  if ($line =~ $searchterm) {$getflag = 1;}
}
