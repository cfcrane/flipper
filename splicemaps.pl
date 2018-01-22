#!/usr/bin/perl
use warnings;
#This script arranges maps of segments of a single linkage group on the basis of a mapped subset of reference markers.
#It is assumed that each marker is unique to a single map in the list of map files.
if (scalar(@ARGV) != 8) {die "This script needs eight arguments.\n";}
open (REF, "< $ARGV[0]"); #miltest_1028a60ABH_002_002_ss3000map.txt
open (MLI, "< $ARGV[1]"); #mapfiles1211.txt
open (OUT, "> $ARGV[2]");
$rscript = $ARGV[3]; 
$rtemp = $ARGV[4];
$spcfgfile = $ARGV[5];
$spoutfile = $ARGV[6];
$keyword = $ARGV[7]; #Final or forced
$getflag = 0;
$k = 0;
while ($line = <REF>) {
  if ($line =~ m/^\n/) {$getflag = 0;}
  if ($getflag == 1) {
    chomp $line;
    #g2276   2276    1580    0.000000        0.010032        0.000000        g206    0.512000        0.689359
    ($marker, $filler) = split(/\s+/, $line, 2);
    $refpos{$marker} = $k;
    $k++;
  }
  if ($line =~ m/$keyword/) {$getflag = 1;}
}
while ($line = <MLI>) {
  chomp $line;
  ($file, $datafile) = split(/\s+/, $line);
  $datafiles{$file} = $datafile;
  open (INP, "< $file");
  open (RSC, "> $rscript");
  $getflag = 0;
  $nlgs = 0;
  $k = 0;
  %localpos = ();
  @localmarkers = ();
  $sumrefposes = 0;
  $nrefposes = 0;
  $lastpos = -1;
  while ($line = <INP>) {
    if ($line =~ m/^\n/) {$getflag = 0;}
    if ($getflag == 1) {
      if ($line =~ m/linkage group/) {$nlgs++;}
      else {
        chomp $line;
        #g85427  2787    21974   0.003000        0.022050        0.002000        g579524 0.548000        0.557427
        @vars = split(/\s+/, $line);
        if (exists($refpos{$vars[0]})) {
          $localpos{$vars[0]} = $k;
          push (@localmarkers, $vars[0]);
          $nrefposes++;
          $sumrefposes += $refpos{$vars[0]};
        }
        if ($lastpos != $vars[4]) {
          push(@{$namesakearrays{$file}}, $vars[0]);
          $nnamesakes{$file}++;
          push(@{$localgroups{$vars[0]}}, $vars[0]);
          $currnamesake = $vars[0];
        }
        else {push(@{$localgroups{$currnamesake}}, $vars[0]);}
        $lastpos = $vars[4];
        $locations{$vars[0]} = $vars[4];
        $maplengths{$file} = $vars[4]; #These files are all ascending, even if marker order is reversed from the reference map.
        $k++;
      }
    }
    if ($line =~ m/$keyword/) {$getflag = 1;}
  }
  if ($nrefposes == 0) {die "nrefposes is 0, something is seriously wrong.\n";}
  $meanrefposes{$file} = $sumrefposes / $nrefposes;
  print RSC "x <- c\($localpos{$localmarkers[0]}";
  for ($i = 1; $i < scalar(@localmarkers); $i++) {print RSC ", $localpos{$localmarkers[$i]}";}
  print RSC "\)\ny <- c\($refpos{$localmarkers[0]}";
  for ($i = 1; $i < scalar(@localmarkers); $i++) {print RSC ", $refpos{$localmarkers[$i]}";}
  print RSC "\)\nsiframe <- data.frame(x, y)\nlmresult = lm(y ~ x, siframe\)\n";
  print RSC "capture.output\(lmresult, file = \"$rtemp\"\)\n";
  close(RSC);
  `R CMD BATCH $rscript`;
  open (RUT, "< $rtemp");
  while ($line = <RUT>) {
    if ($line =~ m/ntercept/) {
      $pine = <RUT>;
      chomp $pine;
      @tars = split(/\s+/, $pine);
      #print OUT "slope = $tars[2] for $file\n";
      $slopes{$file} = $tars[2];
    }
  }
  print OUT "file = $file meanrefposes = $meanrefposes{$file} slopes = $slopes{$file} namesake count = $nnamesakes{$file}\n";
  close (INP);
}
for $key (sort {$meanrefposes{$a} <=> $meanrefposes{$b}} (keys(%meanrefposes))) {
  push(@fileorder, $key);
  if ($slopes{$key} < 0) {
    @temp = ();
    for $i (0 .. $#{$namesakearrays{$key}}) {push(@temp, $namesakearrays{$key}[$i]);}
    @temp = reverse(@temp);
    $namesakearrays{$key} = ();
    push(@{$namesakearrays{$key}}, @temp);
    for ($i = 0; $i < scalar(@temp); $i++) {
      for $j (0 .. $#{$localgroups{$temp[$i]}}) {
        $locations{$localgroups{$temp[$i]}[$j]} = $maplengths{$key} - $locations{$localgroups{$temp[$i]}[$j]};
      }
    }
  }
}
$lastend = 0;
for ($i = 1; $i < scalar(@fileorder); $i++) {
  #print OUT "$fileorder[$i-1] $fileorder[$i] $meanrefposes{$fileorder[$i-1]} $meanrefposes{$fileorder[$i]} $slopes{$fileorder[$i-1]} $slopes{$fileorder[$i]}\n";
  @firstarray = ();
  @secondarray = ();
  @splicearray = ();
  for $j (0 .. $#{$namesakearrays{$fileorder[$i-1]}}) {$firstarray[$j] = $namesakearrays{$fileorder[$i-1]}[$j];}
  for $j (0 .. $#{$namesakearrays{$fileorder[$i]}}) {$secondarray[$j] = $namesakearrays{$fileorder[$i]}[$j];}
  for ($j = -6; $j < 0; $j++) {push(@splicearray, $firstarray[$j]);}
  for ($j = 0; $j < 6; $j++) {push(@splicearray, $secondarray[$j]);}
  #print OUT "splice zone = $splicearray[0]";
  #for ($j = 1; $j < scalar(@splicearray); $j++) {print OUT " $splicearray[$j]";}
  #print OUT "\n";
  open (CFG, "> $spcfgfile");
  print CFG "left_abutment = $splicearray[0]\n";
  for ($j = 1; $j < 11; $j++) {print CFG "marker = $splicearray[$j]\n";}
  print CFG "right_abutment = $splicearray[11]\n";
  print CFG "outfile = $spoutfile\ninputformat = mapmaker\n";
  print CFG "firstdatafile = $datafiles{$fileorder[$i-1]}\nseconddatafile = $datafiles{$fileorder[$i]}\n";
  close(CFG);
  #Call locally permuting C program here.
  `./permutesplice.exe $spcfgfile`;
  $spflag = 0;
  if (-e $spoutfile) {
    @spnames = ();
    @offsets = ();
    open (SPO, "< $spoutfile");
    while ($vine = <SPO>) {
      if ($vine =~ m/^\n/) {$spflag = 0;}
      if ($spflag == 1) {
        #g1000 0.002012
        chomp $vine;
        ($spname, $offset) = split(/\s+/, $vine);
        push(@spnames, $spname);
        push(@offsets, $offset);
      }
      if ($vine =~ m/Verdict/) {$spflag = 1;}
    }
    #`rm $spoutfile`; #Prevent re-use of old splice files at other junctions.
  }
  else {die "Annealed join failed for $fileorder[$i-1] and $fileorder[$i]\n";}
  #print OUT "Order and location of namesake markers:\n";
  #for ($j = 0; $j < scalar(@spnames); $j++) {print OUT "$spnames[$j] $offsets[$j]\n";}
  #print OUT "End of splicing markers.\nBeginning or resuming map:\n";
  if ($i == 1) { #Print out the first map up to left_abutment.
    for ($j = 0; $j < scalar(@firstarray) - 6; $j++) {
      for $k (0 .. $#{$localgroups{$firstarray[$j]}}) {
        print OUT "$localgroups{$firstarray[$j]}[$k]\t$locations{$localgroups{$firstarray[$j]}[$k]}\n";
      }
    }
    $currpos = $locations{$firstarray[-6]};
  }
  else {
    #print OUT "Value of v upon entering i-loop is $v\n";
    $currpos = $v + $locations{$firstarray[-6]} - $locations{$firstarray[-7]};
  }
  for ($j = 0; $j < 12; $j++) {
    $v = $currpos + $offsets[$j];
    for $k (0 .. $#{$localgroups{$spnames[$j]}}) {print OUT "$localgroups{$spnames[$j]}[$k]\t$v\n";}
  }
  $lastend = $v;
  for ($j = 6; $j < scalar(@secondarray) - 6; $j++) {
    for $k (0 .. $#{$localgroups{$secondarray[$j]}}) {
      $v = $lastend + $locations{$localgroups{$secondarray[$j]}[$k]} - $locations{$secondarray[5]};
      print OUT "$localgroups{$secondarray[$j]}[$k]\t$v\n";
    }
  }
  if ($i == scalar(@fileorder) -1) {
    for ($j = scalar(@secondarray) - 6; $j < scalar(@secondarray); $j++) { #Tack on the end of the last map.
      for $k (0 .. $#{$localgroups{$secondarray[$j]}}) {
        $v = $lastend + $locations{$localgroups{$secondarray[$j]}[$k]} - $locations{$secondarray[5]};
        print OUT "$localgroups{$secondarray[$j]}[$k]\t$v\n";
      }
    }
  }
  $lastend = $v;
}
