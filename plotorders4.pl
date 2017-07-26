#!/usr/bin/perl
use warnings;
#This script plots markers positions of one parent versus the other parent.
open (MAP, "< $ARGV[0]"); #pagxxivtestafakfilemap.txt
open (OUT, "> $ARGV[1]"); #pagxxivtestaplotout.txt
open (RSC, "> $ARGV[2]"); #pagxxivtestaplot.R
$pdffilestem = $ARGV[3]; #pagxxivtestaplot
if ($ARGV[4] ne "NULL") {open (KLL, "< $ARGV[4]");} #killedlociintesta.txt
$pdfheight = $ARGV[5]; #4
$pdfwidth = $ARGV[6]; #4
$cutoffcount = $ARGV[7]; #30
if (scalar(@ARGV) != 8) {die "This script needs eight arguments.\n";}
$getflag = 0;
$sensflag = 0;
$timescalled = 0;
%onmap = ();
@mapcounts = ();
while ($line = <MAP>) {
  if ($line =~ m/^\n/) {$getflag = 0; $sensflag = 0;}
  chomp $line;
  if ($line =~ m/removing misfit/) { #This indicates a map stripped of misfit loci using continuation = 1.
    @mapcounts = (); #Clean out values from the first map.
    @claimedcounts = ();
    @mapgroups = ();
    @map = ();
    %killedloci = ();
    %onmap = ();
    %zerohash = ();
    %mapplaces = ();
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
    %killedloci = ();
    %onmap = ();
    %zerohash = ();
    %mapplaces = ();
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
      $basenames{$vars[0]} = $part1;
      if ($part2 == 1) {$sisters{$vars[0]} = $basenames{$vars[0]}."-2";}
      elsif ($part2 == 2) {$sisters{$vars[0]} = $basenames{$vars[0]}."-1";}
      $onmap{$vars[0]} = 1;
      push(@{$map[$lg]}, $vars[0]);
      $mapgroups{$vars[0]} = $lg;
      $place++;
      $mapplaces{$vars[0]} = $place;
      $mapcounts[$lg]++;
      $killedloci{$vars[0]} = 0; #Overwrite killed entries below.
      if ($previd > -0.5) {
        $zerohash{$previd} = 1;
        $zerohash{$vars[0]} = 1;
      }
      if ($vars[2] > -0.5 && $vars[3] == 0.000000) { #Color both partners in a zero-recombining bin.
        if ($vars[0] =~ m/-/) {($previd, $filler) = split(/-/, $vars[0]);}
        else {$previd = $vars[0];}
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
if ($ARGV[4] ne "NULL") {
  while ($line = <KLL>) {
    chomp $line;
    @vars = split/\s+/, $line;
    $name = $vars[-1]."-1";
    $killedloci{$name} = 1;
    $name = $vars[-1]."-2";
    $killedloci{$name} = 1;
  }
}
#This presumes that the Flipper output is ordered descendingly by length of linkage group.
for ($i = 0; $i < scalar(@mapcounts); $i++) {
  if ($mapcounts[$i] < $cutoffcount) {next;}
  $maxtotalsame = 0;
  $savedj = -1;
  for ($j = 0; $j < scalar(@mapcounts); $j++) {
    if ($mapcounts[$j] < $cutoffcount) {next;}
    if ($j == $i) {next;}
    $totalsame = 0;
    for $m (0 .. $#{$map[$i]}) {
      for $n (0 .. $#{$map[$j]}) {
        if ($basenames{$map[$j][$n]} == $basenames{$map[$i][$m]}) {
          $totalsame++;
        }
      }
    }
    if ($totalsame > $maxtotalsame) {$savedj = $j; $maxtotalsame = $totalsame;}
  }
  $pairedwith[$i] = $savedj;
  $pairedwith[$savedj] = $i;
}
for ($i = 0; $i < scalar(@pairedwith); $i++) {
  print OUT "Group $i goes with group $pairedwith[$i] and group $pairedwith[$i] goes with group $pairedwith[$pairedwith[$i]]\n";
}
for ($k = 0; $k < scalar(@pairedwith); $k++) {
  if ($mapcounts[$k] < $cutoffcount) {next;}
  $maxy = 0;
  $ystring = "";
  $pdffile = $pdffilestem.$k.".pdf";
  print RSC "pdf\(\"$pdffile\", height = $pdfheight, width = $pdfwidth\)\nx <- c\(";
  $killedstringx = ""; $mainstringx = ""; $killedstringy = ""; $mainstringy = ""; $zerostringx = ""; $zerostringy = "";
  for ($m = 0; $m < scalar(@pairedwith); $m++) {
    if ($pairedwith[$m] == $k) {
      if ($ystring eq "") {$ystring = $m;}
      else {$ystring .= ", $m";}
      for $i (0 .. $#{$map[$m]}) {
        if ($mapplaces{$map[$m][$i]} > $maxy) {$maxy = $mapplaces{$map[$m][$i]};}
      }
    }
  }
  $sumx = 0; $sumy = 0;
  for $i (0 .. $#{$map[$k]}) {
    if ($killedloci{$map[$k][$i]} == 1 && exists($mapgroups{$sisters{$map[$k][$i]}})) {
      if ($mapgroups{$sisters{$map[$k][$i]}} == $pairedwith[$k]) {
        print OUT "map[$k][$i] = $map[$k][$i] map group = $mapgroups{$map[$k][$i]} sister locus = $sisters{$map[$k][$i]} map group = $mapgroups{$sisters{$map[$k][$i]}}\n";
        if ($killedstringx eq "") {$killedstringx = "$i"; $killedstringy = "$mapplaces{$sisters{$map[$k][$i]}}";}
        else {$killedstringx .= ", $i"; $killedstringy .= ", $mapplaces{$sisters{$map[$k][$i]}}";}
      }
    }
    elsif (exists($zerohash{$map[$k][$i]}) || exists($zerohash{$sisters{$map[$k][$i]}})) {
      if ($zerostringx eq "") {$zerostringx = "$i"; $zerostringy = "$mapplaces{$sisters{$map[$k][$i]}}";}
      elsif (exists($mapplaces{$sisters{$map[$k][$i]}})) {$zerostringx .= ", $i"; $zerostringy .= ", $mapplaces{$sisters{$map[$k][$i]}}";}
    }
    elsif (exists($mapplaces{$sisters{$map[$k][$i]}})) {
      $sumx++; $sumy++;
      if ($mainstringx eq "") {$mainstringx = "$i"; $mainstringy = "$mapplaces{$sisters{$map[$k][$i]}}";}
      else {$mainstringx .= ", $i"; $mainstringy .= ", $mapplaces{$sisters{$map[$k][$i]}}";}
    }
  }
  print OUT "k = $k max y = $maxy\nkilled string x = $killedstringx\nkilled string y = $killedstringy\n";
  print OUT "k = $k zerostring x = $zerostringx\nzerostring y = $zerostringy\n";
  print OUT "k = $k mainstring x = $mainstringx\nmain string y = $mainstringy\n";
  print OUT "sum x = $sumx sum y = $sumy\n";
  print RSC "$mainstringx\)\nxk <- c\($killedstringx\)\ny <- c\($mainstringy\)\nyk <- c\($killedstringy\)\n";
  print RSC "xz <- c\($zerostringx\)\nyz <- c\($zerostringy\)\n";
  $maxy += 10;
  if ($ystring =~ m/,/) {$ystring = "s ".$ystring;}
  else {$ystring = " ".$ystring;}
  if ($ystring =~ m/,/) {
    print RSC "plot\(x, y, ylim = c\(0, $maxy\), type = \"n\", xaxt = \"n\", yaxt = \"n\", xlab = \"Mapped Position in Linkage Group $k\", ylab = \"Mapped Place in Group$ystring\"\)\n";
  }
  else {
    print RSC "plot\(x, y, ylim = c\(0, $maxy\), type = \"n\", xaxt = \"n\", yaxt = \"n\", xlab = \"Mapped Position in Linkage Group $k\", ylab = \"Mapped Position in Linkage Group$ystring\"\)\n";
  }
  print RSC "axis\(side = 1, at = c\(0";
  for ($i = 10; $i < $mapcounts[$k] + 10; $i += 10) {print RSC ", $i";}
  print RSC "\)\)\n";
  print RSC "axis\(side = 2, at = c\(0";
  for ($i = 10; $i < $maxy + 10; $i += 10) {print RSC ", $i";}
  print RSC "\)\)\n";
  print RSC "points\(x, y, pch = 22, col = \"black\", cex = 0.15\)\n";
  print RSC "points\(xk, yk, pch = 22, col = \"red\", cex = 0.15\)\n";
  print RSC "points\(xz, yz, pch = 22, col = \"#33AAFF\", cex = 0.15\)\n";
  print RSC "dev.off\(\)\n";
}
