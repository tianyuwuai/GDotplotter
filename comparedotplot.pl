
print "\n\n Usage :  perl this_script querychr sbjctchr querygff1 objctgff2 pos_or_order queryinitial sbjctinitial filetype1 file1 * filetype2 file2 *    \n\n";

use strict;
use GD;
use GD::Text::Align;

my @querychr = split(/_/, $ARGV[0]);
my @sbjctchr = split(/_/, $ARGV[1]);

my $queryspecies = $ARGV[5];
my $sbjctspecies = $ARGV[6];

my $figure_file = $ARGV[5].$ARGV[0].".v.".$ARGV[6].$ARGV[1].".compare.png";


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
if($ARGV[4] eq "order")
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
elsif($ARGV[4] eq "pos")
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
$align->set_text($ARGV[5]);
$align->draw($frame_width/2 + $left_curb,65,0);
$align->set_text($ARGV[6]);
$align->draw(65,$frame_height/2 + $top_curb,1.57);

print "accumulated chro length is $accumulated_length\n";
my $countnum = 9;
$countnum = drawdotplot($ARGV[7],$ARGV[8],$countnum,$blue);
$countnum = drawdotplot($ARGV[$countnum],$ARGV[$countnum+1],$countnum+2,$red);
print $countnum;
binmode FIG;
print FIG $img -> png;
close($figure_file);


sub drawdotplot()
{
   my($filetype,$filename,$tempnumber,$color) = @_;
   my $tempposx = 0;
   my $tempposy = 0;
   #print $filetype."+++++".$filename."+++++".$tempnumber."\n";
   if($filetype eq "Blast")
   {
   open(IN, $filename) or die "cannot open infile due to $!.\n";
###### parameters for significant hits, for filtering repeats
   my $EVALUE = $ARGV[$tempnumber];
   my $SCORE  = $ARGV[$tempnumber+1];
   my $HITNUM = $ARGV[$tempnumber+2];
   $tempnumber = $tempnumber+3;
   my $isoutput = 0;
   my $lastquery = "";
   my $hitnum = 0;
   while(<IN>)
   {
#print $_;
   $_ =~ s/[\n\r]//g;

   my ($querygene, $sbjctgene, $identity, $matchlen, $mismatchnum, $gaplen, $querystart, $queryend, $sbjctstart, $sbjctend, $evalue, $score) = split(/\t/, $_   );

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

# max hit can not more than HITNUM;
   if($lastquery ne $querygene)
   {
      $hitnum = 1;
      $lastquery = $querygene;
   }
   else
   {
      $hitnum ++;
   }
   if($hitnum > $HITNUM){next;}


# blast score should be more than SCORE
   if($score < $SCORE){next;}

# e-value should be less than EVALUE
   if($evalue > $EVALUE){next;}

   ($tempposx, $tempposy) = getpos($querygene, $sbjctgene, $querychr, $sbjctchr);
   #print $tempposx."*******".$tempposy."\n";
   $img -> filledRectangle($tempposx, $tempposy, $tempposx+1, $tempposy+1, $color);
   }
   }
   elsif($filetype eq "Mcscan")
   {
      open(IN, $filename) or die "cannot open infile due to $!.\n";

###### parameters for significant hits, for filtering repeats
      my $EVALUE = $ARGV[$tempnumber+1];
      my $SCORE  = $ARGV[$tempnumber];
      my $LENGTH = $ARGV[$tempnumber+2];
      $tempnumber = $tempnumber+3;
      my $isoutput = 0;
      my $lastquery = "";      
      my $flag = 0;
      while(<IN>)
      {
      $_ =~ s/[\n\r]//g;

      if($_ eq ""){next;}
   
      if($_ =~ /^## Alignment /){
         $flag = 0;
         my @array = split(/\s+/, $_);
         my @temp = split(/=/, $array[3]);
         my $score = $temp[1];
         @temp = split(/=/, $array[4]);
         my $evalue = $temp[1];
         @temp = split(/=/, $array[5]);
         my $length = $temp[1];
         if($score < $SCORE || $evalue > $EVALUE || $length < $LENGTH){$flag = 1;}
      #print $score." < ".$SCORE."||".$evalue." > ".$EVALUE." || ".$length." < ".$LENGTH."*******".$flag."\n";
      # else {$flag = 0;}
         next;}
   
      if($flag eq 1){next;}

      if($_ =~ /^#/){next;}

      my ($alignnum, $querygene, $sbjctgene, $evalue) = split(/\t/, $_);

      my $querychr = $querygene2chr{$querygene};
      my $sbjctchr = $sbjctgene2chr{$sbjctgene};

# only selected chromosome ;
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
      ($tempposx, $tempposy) = getpos($querygene, $sbjctgene, $querychr, $sbjctchr);
      $img -> filledRectangle($tempposx, $tempposy, $tempposx+1, $tempposy+1, $color);
      #drawdot($querygene, $sbjctgene, $querychr, $sbjctchr,$color);
      }
   }
   elsif($filetype eq "Colinearscan")
   {
      open(IN, $filename) or die "cannot open infile due to $!.\n";
      my $isoutput = 0;
      my $LENGTH = $ARGV[$tempnumber];
      my $PVALUE = $ARGV[$tempnumber+1];
      $tempnumber = $tempnumber+2;
      my $lastquery = "";
      my $number = 0;
      my @AllgeneInABlock = ();
      while(<IN>)
      {
#print $_;
         $_ =~ s/[\n\r]//g;
      	 my @a = split(/\s/, $_);
      	 if($a[0] eq ""){next;}
      	 if($a[0] =~ /^\+/){next;}
      	 if($a[0] eq "the")
      	 {
	     $isoutput = 1;
       	     if($a[4] < $LENGTH){$isoutput = 0; next;}
         }
	 if($isoutput eq 0){next;}
	 if($a[0] =~ /^>LOCAL/)
	 {
	    if($a[3] < $PVALUE)
      	    {
	       for(my $i=0; $i<=$#AllgeneInABlock; $i++)
      	       {
	          my @temparr = split(/\t/, $AllgeneInABlock[$i]);
      		  my $tempx = $temparr[0];
      		  my $tempy = $temparr[1];
		  #my $tempposx, $tempposy = getpos($querygene, $sbjctgene, $querychr, $sbjctchr);
		  #$img -> filledRectangle($tempposx, $tempposy, $tempposx+1, $tempposy+1, $color);
		  $img -> filledRectangle($tempx, $tempy, $tempx+1, $tempy+1, $color); 
               }
            }      
         @AllgeneInABlock = ();
   	 $number = 0;
   	 next;
         }
      	 my $id1 = $a[0];
      	 my $id2 = $a[2];
      	 my $querychr = $querygene2chr{$id1};
      	 my $sbjctchr = $sbjctgene2chr{$id2};

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
         ($tempposx, $tempposy) = getpos($id1, $id2, $querychr, $sbjctchr);
	 $AllgeneInABlock[$number] = $tempposx."\t".$tempposy; 
      	 $number ++;
      }
   }
   elsif($filetype eq "Adhore")
   {
      open(IN, $filename) or die "cannot open infile due to $!.\n";
      my $isoutput = 0;
      my $lastquery = "";
      my $hitnum = 0;
      while(<IN>)
      {
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
	 ($tempposx, $tempposy) = getpos($querygene, $sbjctgene, $querychr, $sbjctchr);
   	 $img -> filledRectangle($tempposx, $tempposy, $tempposx+1, $tempposy+1, $color);
      }
   }
   return $tempnumber;
}
sub getpos()
{
	
#print "query $querychr sbjct $sbjctchr\n";

   my ($posx1, $posy1, $posx2, $posy2, $selfhit1x, $selfhit1y, $selfhit2x, $selfhit2y);
   my ($querygene, $sbjctgene, $querychr, $sbjctchr) = @_;
   if($ARGV[4] eq "order")
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
   elsif($ARGV[4] eq "pos")
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
      }
   }
   #print $posx1."++++".$posy1."\n";
   return ($posx1,$posy1);  

   # $img -> filledRectangle($posx1, $posy1, $posx1+1, $posy1+1, $color);
   
}


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
