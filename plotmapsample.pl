#!/usr/bin/perl
use warnings;
use Math::Random::MT;
open (MAP, "< $ARGV[0]"); #splicemapsout1211.txt
open (PHY, "< $ARGV[1]"); #miltest_1028a60outa.txt
open (OUT, "> $ARGV[2]");
$rscript = $ARGV[3];
open (RSC, "> $rscript");
$pdffile = $ARGV[4];
$iseed = $ARGV[5]; #11279929
$N = $ARGV[6]; #600
$height = $ARGV[7]; #7.5
$width = $ARGV[8]; #7
CORE::srand($iseed);
for ($i = 0; $i < 624; $i++) { #Set up seeds for the Mersenne Twister.  Each seed is an unsigned integer.
  $v = int(CORE::rand(429496));
  for ($j = 0; $j < 4; $j++) {$v .= int(CORE::rand(10));}
  $av[$i] = $v;
}
$rng = Math::Random::MT->new(@av);
$linecount = 0;
while ($line = <MAP>) {
  if ($line !~ m/^file =/) {$linecount++;}
}
$limit = $N / $linecount;
seek(MAP, 0, 0);
while ($line = <MAP>) {
  if ($line !~ m/^file =/) {
    chomp $line;
    @vars = split(/\t/, $line);
    $r = $rng->rand();
    if ($r < $limit) {push(@mappedorder, $vars[0]); $inplot{$vars[0]} = 1;}
  }
}
$telocount = 0;
while ($line = <PHY>) {
  #0 1000000 0.000000 TELO
  #768885 768885 0.000000 1
  chomp $line;
  if ($line =~ m/TELO/) {
    $telocount++;
    if ($telocount == 2) {last;}
  }
  else {
    @vars = split(/\s+/, $line);
    $name = "g".$vars[0];
    if (exists($inplot{$name})) {push(@physorder, $name);}
  }
}
for ($i = 0; $i < scalar(@physorder); $i++) {print OUT "$physorder[$i] $mappedorder[$i]\n";}
print RSC "x <- c\(0";
for ($i = 1; $i < scalar(@physorder); $i++) {print RSC ", $i";}
print RSC "\)\n";
print RSC "y <- c\(";
for ($j = 0; $j < scalar(@mappedorder); $j++) {
  if ($mappedorder[$j] eq $physorder[0]) {print RSC $j;}
}
for ($i = 1; $i < scalar(@physorder); $i++) {
  for ($j = 0; $j < scalar(@mappedorder); $j++) {
    if ($mappedorder[$j] eq $physorder[$i]) {print RSC ", $j";}
  }
}
print RSC "\)\n";
print RSC "pdf\(\"$pdffile\", height = $height, width = $width\)\n";
$ylim = scalar(@mappedorder);
print RSC "plot\(x, y, ylim = c\(0, $ylim\), type = \"n\", xaxt = \"n\", yaxt = \"n\", xlab = \"Physical Marker Order\", ylab = \"Mapped Marker Order\"\)\n";
print RSC "axis\(side = 1, at = c\(0";
for ($i = 50; $i < scalar(@mappedorder); $i += 50) {print RSC ", $i";}
print RSC "\)\)\n";
print RSC "axis\(side = 2, at = c\(0";
for ($i = 50; $i < scalar(@physorder); $i += 50) {print RSC ", $i";}
print RSC "\)\)\n";
print RSC "points\(x, y, pch = 22, col = \"black\", cex = 0.25\)\n";
print RSC "dev.off\(\)\n";
close(RSC);
`/Library/Frameworks/R.framework/Versions/Current/Resources/bin/R CMD BATCH $rscript`;
