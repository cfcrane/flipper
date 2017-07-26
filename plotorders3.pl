#!/usr/bin/perl
use warnings;
#This script plots a separate-parental map order versus the actual locus order for a simulated data set.
open (REF, "< $ARGV[0]"); #pagxxivtestaouta.txt
open (MAP, "< $ARGV[1]"); #pagxxivtestafakfilemap.txt
if ($ARGV[2] ne "NULL") {open (KLL, "< $ARGV[2]");} #killedlociintesta.txt
open (OUT, "> $ARGV[3]"); #pagxxivtestaplotout.txt
open (RSC, "> $ARGV[4]"); #pagxxivtestaplot.R
$pdffilestem = $ARGV[5]; #pagxxivtestaplot
$pdfheight = $ARGV[6]; #4
$pdfwidth = $ARGV[7]; #4
if (scalar(@ARGV) != 8) {die "This script needs eight arguments.\n";}
#We need to go back and forth several times between the same files.
$maxlocus = -1;
while ($line = <REF>) {
  if ($line !~ m/TELO/) {
    chomp $line;
    @vars = split(/\s+/, $line);
    #4012 10009 0.832775 1
    #Below, use key for nnnn-1 and synonyms{key} for nnnn-2.
    if ($vars[1] > $vars[0]) {$synonyms{$vars[0]} = $vars[1];}
    if ($vars[1] > $maxlocus) {$maxlocus = $vars[1];}
  }
}
seek(REF, 0, 0);
if ($ARGV[2] ne "NULL") {
  while ($line = <KLL>) {
    chomp $line;
    @vars = split/\s+/, $line;
    $killedloci{$vars[-1]} = 1;
    $killedloci{$synonyms{$vars[-1]}} = 1;
  }
}
else {%killedloci = ();}
for ($i = 0; $i <= $maxlocus; $i++) {
  $actgroups[$i] = -1;
  $mapgroups[$i] = -1;
  $actplaces[$i] = -1;
  $mapplaces[$i] = -1;
}
$getflag = 0;
$sensflag = 0;
$timescalled = 0;
%onmap = ();
@mapcounts = ();
@claimedcounts = ();
while ($line = <MAP>) {
  if ($line =~ m/^\n/) {$getflag = 0; $sensflag = 0;}
  chomp $line;
  if ($line =~ m/removing misfit/) { #This indicates a map stripped of misfit loci using continuation = 1.
    @mapcounts = (); #Clean out values from the first map.
    @claimedcounts = ();
    @mapgroups = ();
    @map = ();
    %onmap = ();
    %zerohash = ();
    @mapplaces = ();
    $sensflag = 1; #Here we go again.
    $getflag = 0;
    $timescalled++;
  }
  if ($line =~ m/forced placement/) { #This indicates a second map produced with continuation = 2.
    #This block will follow the continuation == 1 block if both are printed.
    @mapcounts = (); #Clean out values from the descending or misfit map.
    @claimedcounts = ();
    @mapgroups = ();
    @map = ();
    %onmap = ();
    %zerohash = ();
    @mapplaces = ();
    $sensflag = 1; #Here we go again.
    $getflag = 0;
    $timescalled++;
  }
  if ($sensflag == 1) {
    if ($line =~ m/linkage group/) {
      $previd = -1;
      $place = -1;
      @vars = split(/\s+/, $line);
      #linkage group 0 start = 830     end     1199    count = 121     length = 1.087500
      $lg = $vars[2];
      $mapcounts[$lg] = 0;
      $claimedcounts[$lg] = $vars[10];
      $getflag = 1;
    }
    elsif ($getflag == 1) {
      @vars = split(/\s+/, $line);
      #combined: 2392    2392    2933    0.010000        0.145000        0.010000        0.545000        1.2308
      #separate: 3236-1  6472    1354    0.010000        0.035075        0.005000        2433-1  0.540000        0.3182
      ($part1, $part2) = split(/-/, $vars[0]);
      if ($vars[0] =~ m/-1$/) {$usedname = $part1;}
      elsif ($vars[0] =~ m/-2$/) {$usedname = $synonyms{$part1};}
      $onmap{$usedname} = 1;
      push(@{$map[$lg]}, $usedname);
      $mapgroups[$usedname] = $lg;
      $place++;
      $mapplaces[$usedname] = $place;
      $mapcounts[$lg]++;
      if ($previd > -0.5) {
        $zerohash{$previd} = 1;
        $zerohash{$usedname} = 1;
      }
      if ($vars[2] > -0.5 && $vars[3] == 0.000000) { #Color both partners in a zero-recombining bin.
        $previd = $usedname;
      }
      else {$previd = -1;}
    }
  }
  if ($line =~ m/descending/) {$sensflag = 1;} #Set sensflag after lines that key on "linkage group", since the "descending" line contains the string "linkage group".
}
print OUT "The so-called final map was seen $timescalled times.\n";
for ($i = 0; $i < scalar(@mapcounts); $i++) {print OUT "mapcounts for group $i = $mapcounts[$i]\n";}
for ($i = 0; $i < scalar(@mapcounts); $i++) {
  if ($mapcounts[$i] != $claimedcounts[$i]) {
    print OUT "Claimed and counted markers in linkage group $lg differ: claimed = $claimedcounts[$i] counted = $mapcounts[$i]\n";
  }
}
$lg = -1;
while ($line = <REF>) {
  chomp $line;
  if ($line =~ m/TELO/) {$lg++; $actcounts[$lg] = 0; $place = -1;}
  else {
    #2692 8689 0.999126 1
    @vars = split(/\s+/, $line);
    if (exists($onmap{$vars[1]})) {
      push(@{$orders[$lg]}, $vars[1]);
      $actgroups[$vars[1]] = $lg;
      $place++;
      $actplaces[$vars[1]] = $place;
    }
    $actcounts[$lg]++;
  }
}
#Diagnostics added 10-1-2016.
for ($i = 0; $i < scalar(@map); $i++) {
  for $j (0 .. $#{$map[$i]}) {
    print OUT "DIAG: linkage group i = $i locus j = $j map[$i][$j] = $map[$i][$j] onmap{$map[$i][$j]} = $onmap{$map[$i][$j]} actgroups[$map[$i][$j]] = $actgroups[$map[$i][$j]] actplaces[$map[$i][$j]] = $actplaces[$map[$i][$j]]\n";
  }
}
for ($i = 0; $i < scalar(@map); $i++) {
  $abssum = 0;
  $altsum = 0;
  for $j (0 .. $#{$map[$i]}) {
    $abssum += abs($actplaces[$map[$i][$j]] - $mapplaces[$map[$i][$j]]);
    $altsum += abs($actplaces[$map[$i][$j]] - ($mapcounts[$i] - $mapplaces[$map[$i][$j]] - 1));
  }
  print OUT "abssum for group $i = $abssum alternate sum = $altsum\n";
  if ($abssum > $altsum) {
    for $j (0 .. $#{$map[$i]}) {$mapplaces[$map[$i][$j]] = $mapcounts[$i] - $mapplaces[$map[$i][$j]] - 1;}
    for $j (0 .. $#{$map[$i]}) {$printmap[$i][$j] = $map[$i][$mapcounts[$i]-$j-1];}
  }
  else {
    for $j (0 .. $#{$map[$i]}) {$printmap[$i][$j] = $map[$i][$j];}
  }
}
print OUT "Actual and mapped locus counts by actual and mapped linkage group:\n";
#The following will give warnings when there are more groups on the recombinational map than in the actual physical map.
for ($i = 0; $i < scalar(@mapcounts); $i++) {print OUT "group $i: $actcounts[$i] and $mapcounts[$i] markers\n";}
print OUT "Actual order:\n";
for ($i = 0; $i < scalar(@orders); $i++) {
  for $j (0 .. $#{$orders[$i]}) {print OUT "$i $j $orders[$i][$j]\n";}
}
print OUT "MAP with mapped and actual group memberships and locus positions:\n";
for ($i = 0; $i < scalar(@map); $i++) {
  for $j (0 .. $#{$map[$i]}) {
    print OUT "$i $j $map[$i][$j] $mapgroups[$map[$i][$j]] $mapplaces[$map[$i][$j]] $actgroups[$map[$i][$j]] $actplaces[$map[$i][$j]]\n";
  }
}
#The following will cause warnings when there are more mapped than actual linkage groups.
for ($k = 0; $k < scalar(@mapcounts); $k++) { #mapcounts has an extra void lg for the last one
  $pdffile = $pdffilestem.$k.".pdf";
  print RSC "pdf\(\"$pdffile\", height = $pdfheight, width = $pdfwidth\)\nx <- c\(";
  $killedstringx = ""; $mainstringx = ""; $killedstringy = ""; $mainstringy = ""; $zerostringx = ""; $zerostringy = "";
  $maxy = 0;
  for ($i = 0; $i < $mapcounts[$k]; $i++) {
    if ($actplaces[$printmap[$k][$i]] > $maxy) {$maxy = $actplaces[$printmap[$k][$i]];}
    if (exists($killedloci{$printmap[$k][$i]})) {
      if ($killedstringx eq "") {$killedstringx = "$i"; $killedstringy = "$actplaces[$printmap[$k][$i]]";}
      else {$killedstringx .= ", $i"; $killedstringy .= ", $actplaces[$printmap[$k][$i]]";}
    }
    elsif (exists($zerohash{$printmap[$k][$i]})) {
      if ($zerostringx eq "") {$zerostringx = "$i"; $zerostringy = "$actplaces[$printmap[$k][$i]]";}
      else {$zerostringx .= ", $i"; $zerostringy .= ", $actplaces[$printmap[$k][$i]]";}
    }
    else {
      if ($mainstringx eq "") {$mainstringx = "$i"; $mainstringy = "$actplaces[$printmap[$k][$i]]";}
      else {$mainstringx .= ", $i"; $mainstringy .= ", $actplaces[$printmap[$k][$i]]";}
    }
  }
  print RSC "$mainstringx\)\nxk <- c\($killedstringx\)\ny <- c\($mainstringy\)\nyk <- c\($killedstringy\)\n";
  print RSC "xz <- c\($zerostringx\)\nyz <- c\($zerostringy\)\n";
  $maxy += 10;
  print RSC "plot\(x, y, ylim = c\(0, $maxy\), type = \"n\", xaxt = \"n\", yaxt = \"n\", xlab = \"Mapped Position in Linkage Group $k\", ylab = \"Actual Position\"\)\n";
  print RSC "axis\(side = 1, at = c\(0";
  for ($i = 10; $i < $mapcounts[$k] + 10; $i += 10) {print RSC ", $i";}
  print RSC "\)\)\n";
  print RSC "axis\(side = 2, at = c\(0";
  for ($i = 10; $i < $actcounts[$k] + 10; $i += 10) {print RSC ", $i";}
  print RSC "\)\)\n";
  print RSC "points\(x, y, pch = 22, col = \"black\", cex = 0.15\)\n";
  print RSC "points\(xk, yk, pch = 22, col = \"red\", cex = 0.15\)\n";
  print RSC "points\(xz, yz, pch = 22, col = \"#33AAFF\", cex = 0.15\)\n";
}
print RSC "dev.off\(\)\n";
