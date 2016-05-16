### input files:
#   blastout between genes
#   chromosome length files: 1 for within genome blast, 2 for between genome blast
#   gene positions in bp:    1 for within genome blast, 2 for between genome blast
###
print "\n\n Usage :  perl this_script querychr sbjctchr querygff1 objctgff2 adhore pos_or_order queryinitial sbjctinitial\n\n";

print "some parameters may need to be changed int this_script\n\n";

use strict;
use GD;
use GD::Text::Align;

my @querychr = split(/_/, $ARGV[0]);
my @sbjctchr = split(/_/, $ARGV[1]);

my $queryspecies = $ARGV[6];
my $sbjctspecies = $ARGV[7];

my $figure_file = $ARGV[6].$ARGV[0].".v.".$ARGV[7].$ARGV[1].".adhore2dotplot.png";
#my $chrolenfile1 = $ARGV[2];
#my %sbjctchr2len;
#input_chro_len1($chrolenfile1);

#foreach my $chr(sort(keys(%sbjctchr2len)))
#{
#   print "$chr $sbjctchr2len{$chr} \n";
#}

###### parameters for significant hits, for filtering repeats
#my $EVALUE = 1e-5;
#my $SCORE  = $ARGV[9];
#my $HITNUM = $ARGV[10];
#my $REPNUM = 20;

###### parameters for input files
my $querygenegff = $ARGV[2]; ### for gne positions
open(QGFF, $querygenegff) or die "cannot open $querygenegff due to $!\n";

my %querygene2pos;
my %querygene2order;
my %querygene2chr;
my %querychr2plen;
my %querychr2olen;
my $geneNumOnAchr;
my $lastchr = "";
while(<QGFF>)
{
  my @a = split("\t", $_);
  if($a[0] eq $lastchr)
  {
    $geneNumOnAchr ++;
  }
  else
  {
    $geneNumOnAchr = 1;
  }
  $querygene2order{$a[1]}=$geneNumOnAchr;
  $querygene2pos{$a[1]} = $a[2];
  my $chr = $a[0];
  $querygene2chr{$a[1]} = $chr;


  if($querychr2olen{$a[0]} !~ /^\d/)
  {
     $querychr2olen{$a[0]} = $geneNumOnAchr;
  }
  elsif($querychr2olen{$a[0]} < $geneNumOnAchr)
  {
     $querychr2olen{$a[0]} = $geneNumOnAchr;
  }

  if($querychr2plen{$a[0]} !~ /^\d/)
  {
     $querychr2plen{$a[0]} = $a[2];
  }
  elsif($querychr2plen{$a[0]} < $a[2])
  {
     $querychr2plen{$a[0]} = $a[2];
  }
  $lastchr = $chr;
  next;
}

my $sbjctgenegff = $ARGV[3]; ### for gne positions
open(SGFF, $sbjctgenegff) or die "cannot open $sbjctgenegff due to $!\n";

my %sbjctgene2pos;
my %sbjctgene2order;
my %sbjctgene2chr;
my %sbjctchr2plen;
my %sbjctchr2olen;
   $geneNumOnAchr = 0;
   $lastchr = "";
while(<SGFF>)
{
  my @a = split("\t", $_);
  if($a[0] eq $lastchr)
  {
    $geneNumOnAchr ++;
  }
  else
  {
    $geneNumOnAchr = 1;
  }
  $sbjctgene2order{$a[1]}=$geneNumOnAchr;
  $sbjctgene2pos{$a[1]} = $a[2];
  my $chr = $a[0];
  $sbjctgene2chr{$a[1]} = $chr;


  if($sbjctchr2olen{$a[0]} !~ /^\d/)
  {
     $sbjctchr2olen{$a[0]} = $geneNumOnAchr;
  }
  elsif($sbjctchr2olen{$a[0]} < $geneNumOnAchr)
  {
     $sbjctchr2olen{$a[0]} = $geneNumOnAchr;
  }

  if($sbjctchr2plen{$a[0]} !~ /^\d/)
  {
     $sbjctchr2plen{$a[0]} = $a[2];
  }
  elsif($sbjctchr2plen{$a[0]} < $a[2])
  {
     $sbjctchr2plen{$a[0]} = $a[2];
  }
  $lastchr = $chr;
}
#foreach my $chr(sort(keys(%querychr2len)))
#{
#  print "$chr $querychr2len{$chr}\n";
#}

###### parameters of the figures
my $cell_width = 2;
my $cell_height= $cell_width;
my $frame_width = 2000;
my $frame_height= 2000;

my $left_curb = 200;
my $top_curb  = 200;

open(FIG, ">".$figure_file) or die "cannot open $figure_file due to $!\n";


my @querychrlen = ();
my %querychr2order;
my @sbjctchrlen = ();
my %sbjctchr2order;
if($ARGV[5] eq "order")
{
   for(my $i=0; $i<=$#querychr; $i++)
   {
      my $chr = $querychr[$i];
      $querychrlen[$#querychrlen+1] = $querychr2olen{$chr};
      $querychr2order{$chr} = $#querychrlen;
   }
   
   for(my $i=0; $i<=$#sbjctchr; $i++)
   {
      my $chr = $sbjctchr[$i];
      $sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2olen{$chr};
      $sbjctchr2order{$chr} = $#sbjctchrlen;
   }
}
elsif($ARGV[5] eq "pos")
{
   for(my $i=0; $i<=$#querychr; $i++)
   {
      my $chr = $querychr[$i];
      $querychrlen[$#querychrlen+1] = $querychr2plen{$chr};
      $querychr2order{$chr} = $#querychrlen;
   }

   for(my $i=0; $i<=$#sbjctchr; $i++)
   {
      my $chr = $sbjctchr[$i];
      $sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2plen{$chr};
      $sbjctchr2order{$chr} = $#sbjctchrlen;
   }
}


my ($genome1_length, $genome2_length) = (0, 0);

$genome1_length = sum_chr_len(@querychrlen);   
$genome2_length = sum_chr_len(@sbjctchrlen);

my $querychrno = $#querychrlen + 1;
my $sbjctchrno = $#sbjctchrlen + 1;

###### calculate scale ratio
my $scale_ratio1 = $genome1_length/$frame_width; ## horizontal
my $scale_ratio2 = $genome2_length/$frame_height;## vertical

###### draw the frame and the saprating lines corresponding to chromosome borders
my $img = GD::Image -> new($frame_width + 3*$left_curb/2, $frame_height + 3*$top_curb/2);
my $white = $img -> colorAllocate(255, 255, 255);
#$img -> transparent($white);
$img -> interlaced('true');

my $black = $img->colorAllocate(0,0,0);
my $red = $img->colorAllocate(255,0,0);
my $blue = $img->colorAllocate(0,0,255);
my $green = $img->colorAllocate(255,0,255);
my $yellow = $img->colorAllocate(250,250,200);
my $gray = $img->colorAllocate(100,100, 100);
print "outerframe: $frame_width, $frame_height\n";

$img -> rectangle($left_curb, $top_curb, $frame_width + $left_curb, $frame_height + $top_curb, $black);
$img -> rectangle($left_curb-1, $top_curb-1, $frame_width + $left_curb+1, $frame_height + $top_curb+1, $black);

my $align = GD::Text::Align->new($img, valign => 'center', halign => 'center', color => $black);
$align->set_font('Arial.ttf',34);

my @query_chro_pos = ();
my @sbjct_chro_pos = ();

my $accumulated_length = 0;
for(my $i=0; $i<=$#querychrlen; $i++)
{
   $accumulated_length += $querychrlen[$i];
   my $length = int($querychrlen[$i]/$scale_ratio1);
   my $posx1 = $left_curb + int($accumulated_length/$scale_ratio1);
   $query_chro_pos[$i] = $posx1;
   my $posy1 = $top_curb;
   my $posx2 = $posx1;
   my $posy2 = $top_curb + $frame_height;

#print "x chromosome: ".$posx1." ".$posy1." ".$posx2." ".$posy2."\n";
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1+1, $posy1, $posx2+1, $posy2, $black);

   my $chr = $querychr[$i];
   $align->set_text($chr);
   $align->draw($posx1-int($length/2), 150,0);
   #$img -> string(gdGiantFont, $posx1-int($length/2)-40, 50, $chr, $red);
}

$accumulated_length = 0;
for(my $i=0; $i<=$#sbjctchrlen; $i++)
{
   $accumulated_length += $sbjctchrlen[$i];
   my $length = int($sbjctchrlen[$i]/$scale_ratio2);
   my $posy1 = $top_curb + int($accumulated_length/$scale_ratio2);
   $sbjct_chro_pos[$i] = $posy1;
   my $posx1 = $left_curb;
   my $posy2 = $posy1;
   my $posx2 = $left_curb + $frame_width;

#print "y chromosome: ".$posx1." ".$posy1." ".$posx2." ".$posy2."\n";
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1, $posy1+1, $posx2, $posy2+1, $black);
   my $chr = $sbjctchr[$i];
   $align->set_text($chr);
   $align->draw(150, $posy1-int($length/2), 1.57);
   #$img -> string(gdGiantFont, 40, $posy1-int($length/2)-40, $chr, $red);
}

$align->set_font('Arial.ttf',70);
$align->set_text($ARGV[6]);
$align->draw($frame_width/2 + $left_curb,65,0);
$align->set_text($ARGV[7]);
$align->draw(65,$frame_height/2 + $top_curb,1.57);

print "accumulated chro length is $accumulated_length\n";

open(IN, $ARGV[4]) or die "cannot open infile due to $!.\n";

my $isoutput = 0;
my $lastquery = "";
my $hitnum = 0;
while(<IN>)
{
#print $_;
   $_ =~ s/[\n\r]//g;
 
   my ($id, $c1, $multiplicon, $querygene, $sbjctgene, $code) = split(/\t/, $_);
   my $tempgene="";
   if(exists $querygene2chr{$sbjctgene} && exists $sbjctgene2chr{$querygene})
   {
     $tempgene=$querygene;
     $querygene=$sbjctgene;
     $sbjctgene=$tempgene;
   }
   my $querychr = $querygene2chr{$querygene};
   my $sbjctchr = $sbjctgene2chr{$sbjctgene};

# only selected chromosome;
   my $is2skip1 = 1;
   for(my $i=0; $i<=$#querychr; $i++)
   {
     if($querychr eq $querychr[$i]){$is2skip1 = 0; last;}
   }
   my $is2skip2 = 1;
   for(my $i=0; $i<=$#sbjctchr; $i++)
   {
     if($sbjctchr eq $sbjctchr[$i]){$is2skip2 = 0; last;}
   }

   if($is2skip1 eq 1 || $is2skip2 eq 1){next;}


#print "query $querychr sbjct $sbjctchr\n";

   my ($posx1, $posy1, $posx2, $posy2, $selfhit1x, $selfhit1y, $selfhit2x, $selfhit2y);

   if($ARGV[5] eq "order")
   {

#print "query to order $querychr2order{$querychr}\n";
      if($querychr2order{$querychr} eq 0)
      {
         $posx1 = $left_curb + $querygene2order{$querygene}/$scale_ratio1;
      }
      else
      {
         $posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2order{$querygene}/$scale_ratio1;
      }

#print "sbjct to order $sbjctchr2order{$sbjctchr}\n";
      if($sbjctchr2order{$sbjctchr} eq 0)
      {
         $posy1 = $top_curb + $sbjctgene2order{$sbjctgene}/$scale_ratio2;
      }
      else
      {
         $posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2order{$sbjctgene}/$scale_ratio2;
      }
   }
   elsif($ARGV[5] eq "pos")
   {

#print "query to order $querychr2order{$querychr}\n";
      if($querychr2order{$querychr} eq 0)
      {
         $posx1 = $left_curb + $querygene2pos{$querygene}/$scale_ratio1;
      }
      else
      {
         $posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2pos{$querygene}/$scale_ratio1;
      }

#print "sbjct to order $sbjctchr2order{$sbjctchr}\n";
      if($sbjctchr2order{$sbjctchr} eq 0)
      {
         $posy1 = $top_curb + $sbjctgene2pos{$sbjctgene}/$scale_ratio2;
      }
      else
      {
         $posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2pos{$sbjctgene}/$scale_ratio2;
      }   }

#print "posx1 $posx1 poxy1 $posy1\n";

   my $color = $gray;


   $img -> filledRectangle($posx1, $posy1, $posx1+1, $posy1+1, $color);
}

binmode FIG;
print FIG $img -> png;

close($figure_file);


sub sum_chr_len()
{
   my @chrlen = @_;
#print "@chrlen[0..$#chrlen]\n";

   my $sumlen = $chrlen[0];
   for(my $i=1; $i<=$#chrlen; $i++)
   {
#    print "input chro len ".$chrlen[$i]."\n";
    $sumlen = $sumlen + $chrlen[$i];
   }
   return $sumlen;
}

