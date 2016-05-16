#########
# files format:
##### gff:
# chr\tgeneid\tstart_pos\tend_pos
# hv1 Hv1g00001 4070  4229
# 4G  Hv1g00002 4311  4498
# 07g Hv1g00003 4596  4787
##### M8:
# the same as BLAST M8 result
##########
# in this script, user can select special chromosome, filter Blast M8 file, set length of each speci-grid-line,
# here, join the chromosome number with "_"
# I also put some parameters ahead the script, for example colors array, input file names. It can be changed conveniently.

# It is a example input command below:
# perl blast2dotplot.pl -i [inputfilename] -o [outputfilename] -fq [querygfffilename] -fs [sbjctgfffilename] -e 1e-5 -s 200 -n 10 -gs 10000000 -gq 10000000 -cs 1_3_6 -cq 2_5 -dp 1 -p n
####### parameters:
# -i  input file
# -fq gff file of query
# -fs gff file of sbjct
# -o  output figure
# -e  e-value of BLAST  default:5e-2
# -s  score of BLAST  default:0
# -n  hitnumber of BLAST  default:50
# -gs grid line length of subject,if none,input "0" default:0
# -gq grid line lenght of query,if none,input "0" default:0
# -cs selected subject chrmosome number default:all
# -cq selected query chromosome number  default:all
# -dp dot's pix in dotplot default:2
# -p  [y|n] draw by position or order default:y
############################################

use strict;
use Getopt::Long;
use GD;
use GD::Text::Align;
# set parameter of figure
my $frame_width=2000;
my $frame_height=2000;
my $left_curb=30;
my $top_curb=30;
my $img=GD::Image -> new($frame_width+1.5*$left_curb,$frame_height+1.5*$top_curb);


########################### total parameters below
##set variate
#files
my $gff_sfile="hv.new.gff"; #gff file
my $gff_qfile="os.new.gff"; #gff file
my $plot_file="os.hv.cds.blastn.m8.new"; #blast file for dotplot
my $output_fig="hv.draw.repeat.jpg"; #output figure
#special parameters
my $queryspe="Q";
my $sbjctspe="S";
#set color
my $white= $img->colorAllocate(255,255,255);
my $black= $img->colorAllocate(0,0,0);
my $gray= $img->colorAllocate(100,100,100);
my $red= $img->colorAllocate(255,0,0);
my $green= $img->colorAllocate(0,255,0);
my $blue= $img->colorAllocate(0,0,255);
my $mintcream= $img->colorAllocate(245,255,250);
my $dodgerblue= $img->colorAllocate(30,144,255);
my $darkviolet= $img->colorAllocate(148,0,211);
my $orange= $img->colorAllocate(255,165,0);
#font
my $align = GD::Text::Align->new($img, valign => 'center', halign => 'center', color => $black);
$align->set_font('Arial.ttf',34);
#screen input
my $_evalue=5e-2; #blast e-value for dotplot
my $_score=0; #blast score for dotplot filter
my $_hitnum=50; #blast hitnumber for dotplot
my $grid_s=0; #sbjct-speci-grid length
my $grid_q=0; #query-speci-grid length
my $chrn_s="all"; #selected chromosome of sbjct
my $chrn_q="all"; #selected chromosome of query
my $dot_pix=2; #dot pix in dotplot
my $poro="y";
#get variate
Getopt::Long::GetOptions(
  'e=f' => \$_evalue,
  's=i' => \$_score,
  'n=i' => \$_hitnum,
  'gs=i' => \$grid_s,
  'gq=i' => \$grid_q,
  'cs=s' => \$chrn_s,
  'cq=s' => \$chrn_q,
  'dp=i' => \$dot_pix,
  'p=s' => \$poro);
########################################

###################### start dotplot part
####### start read gff file
## read query gff and get gene position, then sort hash and get gene order
open(QGFF,$gff_qfile) or die "can't open query gff file.\n";
my %querygene2pos; #gene id -> pos
my %querygene2order; #gene id -> order
my %querygene2chr; #gene id -> chr number
my %querychr2plen; #chr -> chr total pos
my %querychr2olen; #chr -> chr total order
my %tempsort; # temp has for sorting to get order
#read query gff and get gene pos, chr
while (<QGFF>){
  my @a=split("\t",$_);
  if($a[2]>$a[3]){my $t=$a[3];$a[3]=$a[2];$a[2]=$t;} #make $a[3] as max
  $querygene2pos{$a[1]} = $a[2]; # gene->pos
  $a[0] =~ s/[^0-9]//g;
  my $chr=int($a[0]);
  $querygene2chr{$a[1]} = $chr; #gene->chr
  $tempsort{$a[1]}=$chr."g".$a[2]; #temp for sorting
  #chr->total
  if($querychr2plen{$chr} !~ /^\d/){$querychr2plen{$chr} = $a[3];}
  elsif($querychr2plen{$chr} < $a[3]){$querychr2plen{$chr} = $a[3];}
}
close($gff_qfile);
#sort hash and get order
my $geneNumOnAchr;
my $lastchr = "";
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort){
  my @a=split('g',$tempsort{$key});
  if($a[0] eq $lastchr){$geneNumOnAchr ++;}
    else{$geneNumOnAchr = 1;}
    $querygene2order{$key}=$geneNumOnAchr; #gene->order
    #chr->order
    if($querychr2olen{$a[0]} !~ /^\d/){$querychr2olen{$a[0]} = $geneNumOnAchr;}
    elsif($querychr2olen{$a[0]} < $geneNumOnAchr){$querychr2olen{$a[0]} = $geneNumOnAchr;}
    $lastchr=$a[0];
}
undef(%tempsort);

## read sbjct gff and get gene position, then sort hash and get gene order
open(SGFF,$gff_sfile) or die "can't open sbjct gff file.\n";
my %sbjctgene2pos; #gene id -> pos
my %sbjctgene2order; #gene id -> order
my %sbjctgene2chr; #gene id -> chr number
my %sbjctchr2plen; #chr -> chr total pos
my %sbjctchr2olen; #chr -> chr total order
my %tempsort; # temp has for sorting to get order
#read query gff and get gene pos, chr
while (<SGFF>){
  my @a=split("\t",$_);
  if($a[2]>$a[3]){my $t=$a[3];$a[3]=$a[2];$a[2]=$t;} #make $a[3] as max
  $sbjctgene2pos{$a[1]} = $a[2]; # gene->pos
  $a[0] =~ s/[^0-9]//g;
  my $chr=int($a[0]);
  $sbjctgene2chr{$a[1]} = $chr; #gene->chr
  $tempsort{$a[1]}=$chr."g".$a[2]; #temp for sorting
  #chr->total
  if($sbjctchr2plen{$chr} !~ /^\d/){$sbjctchr2plen{$chr} = $a[3];}
  elsif($sbjctchr2plen{$chr} < $a[3]){$sbjctchr2plen{$chr} = $a[3];}
}
close($gff_sfile);
#sort hash and get order
my $geneNumOnAchr;
my $lastchr = "";
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort){
  my @a=split('g',$tempsort{$key});
  if($a[0] eq $lastchr){$geneNumOnAchr ++;}
    else{$geneNumOnAchr = 1;}
    $sbjctgene2order{$key}=$geneNumOnAchr; #gene->order
    #chr->order
    if($sbjctchr2olen{$a[0]} !~ /^\d/){$sbjctchr2olen{$a[0]} = $geneNumOnAchr;}
    elsif($sbjctchr2olen{$a[0]} < $geneNumOnAchr){$sbjctchr2olen{$a[0]} = $geneNumOnAchr;}
    $lastchr=$a[0];
}
undef(%tempsort);
####### end read gff file


############### get selected chromosome total length and order
my @querychrlen = ();
my %querychr2order;
#if selected chr is all chr, change chrn_ to all chr number
if ($chrn_q eq "all"){
  my @a=();
  foreach my $key (sort{$a<=>$b} keys %querychr2plen){push(@a,$key);}
  $chrn_q=join("_",@a);
}
my @querychr=split('_',$chrn_q);

my @sbjctchrlen = ();
my %sbjctchr2order;
#if selected chr is all chr, change chrn_ to all chr number
if ($chrn_s eq "all"){
  my @a=();
  foreach my $key (sort{$a<=>$b} keys %sbjctchr2plen){push(@a,$key);}
  $chrn_s=join("_",@a);
}
my @sbjctchr=split('_',$chrn_s);
if($poro eq "y"){
  for(my $i=0; $i<=$#querychr; $i++){
    $querychr[$i]=~s/[^0-9]//g;
    my $chr = int($querychr[$i]);
    $querychrlen[$#querychrlen+1] = $querychr2plen{$chr}; #sum length
    $querychr2order{$chr} = $#querychrlen; #total order
  }

  for(my $i=0; $i<=$#sbjctchr; $i++){
    $sbjctchr[$i]=~s/[^0-9]//g;
    my $chr = int($sbjctchr[$i]);
    $sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2plen{$chr}; #sum length
    $sbjctchr2order{$chr} = $#sbjctchrlen; #total order
  }
}
else{
  for(my $i=0; $i<=$#querychr; $i++){
    $querychr[$i]=~s/[^0-9]//g;
    my $chr = int($querychr[$i]);
    $querychrlen[$#querychrlen+1] = $querychr2olen{$chr}; #sum length
    $querychr2order{$chr} = $#querychrlen; #total order
  }

  for(my $i=0; $i<=$#sbjctchr; $i++){
    $sbjctchr[$i]=~s/[^0-9]//g;
    my $chr = int($sbjctchr[$i]);
    $sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2olen{$chr}; #sum length
    $sbjctchr2order{$chr} = $#sbjctchrlen; #total order
  }
}
## get sum length
my ($genome1_length, $genome2_length) = (sum_chr_len(@querychrlen), sum_chr_len(@sbjctchrlen));
## get total selected chromosome quantity
my $querychrno = $#querychrlen + 1;
my $sbjctchrno = $#sbjctchrlen + 1;

###### calculate scale ratio
my $scale_ratio1 = $genome1_length/($frame_width); ## horizontal
my $scale_ratio2 = $genome2_length/($frame_height);## vertical

###############################OK, finished all parameters setting

############################### Let's draw dotplot
# draw the frame of dotplot part and the saprating lines corresponding to chromosome borders
$img -> interlaced('true');
#frame
$img -> rectangle($left_curb*1.25,$top_curb*1.25,$left_curb*1.25+$frame_width,$top_curb*1.25+$frame_height,$black);
$img -> rectangle($left_curb*1.25-1,$top_curb*1.25-1,$left_curb*1.25+$frame_width+1,$top_curb*1.25+$frame_height+1,$black);
#lines
my @query_chro_pos = ();
my @sbjct_chro_pos = ();

my $accumulated_length = 0;
for(my $i=0; $i<=$#querychrlen; $i++){
   my $posy1 = $top_curb*1.25;   
   my $posy2 = $top_curb*1.25 + $frame_height;
#draw grid line
   if($grid_q != 0){     
     for(my $sum_len1 = $grid_q;$sum_len1<$querychrlen[$i];$sum_len1+=$grid_q){
       my $posx0 = $left_curb*1.25 + int(($accumulated_length+$sum_len1)/$scale_ratio1);
       $img -> line($posx0, $posy1, $posx0, $posy2, $green);
     }
   }
   #draw lines between chromosome
   $accumulated_length += $querychrlen[$i];
   my $length = int($querychrlen[$i]/$scale_ratio1);
   my $posx1 = $left_curb*1.25 + int($accumulated_length/$scale_ratio1);
   $query_chro_pos[$i] = $posx1;
   my $posx2 = $posx1;
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1+1, $posy1, $posx2+1, $posy2, $black);
   # plot chromosome number
   my $chr = $querychr[$i];
   $align->set_text($queryspe.$chr);
   $align->draw($posx1-int($length/2),50,0);   
}
#sbjct lines and grid lines
$accumulated_length = 0;
for(my $i=0; $i<=$#sbjctchrlen; $i++){
   my $posx1 = $left_curb*1.25;
   my $posx2 = $left_curb*1.25+$frame_width;
   #grid lines
   if($grid_s != 0){     
     for(my $sum_len2 = $grid_s;$sum_len2<$sbjctchrlen[$i];$sum_len2+=$grid_s){
       my $posy0 = $top_curb*1.25 + int(($accumulated_length+$sum_len2)/$scale_ratio2);
       $img -> line($posx1, $posy0, $posx2, $posy0, $green);
     }
   }
   #lines between chromosomes
   $accumulated_length += $sbjctchrlen[$i];
   my $length = int($sbjctchrlen[$i]/$scale_ratio2);
   my $posy1 = $top_curb*1.25 + int($accumulated_length/$scale_ratio2);
   $sbjct_chro_pos[$i] = $posy1;
   my $posy2 = $posy1;
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1, $posy1+1, $posx2, $posy2+1, $black);
   #chromosome numbers
   my $chr = $sbjctchr[$i];
   $align->set_text($sbjctspe.$chr);
   $align->draw(50, $posy1-int($length/2), 1.57);
}
###########
###########start draw dot
open(DOT,$plot_file) or die "can't open blast file for dotplot.\n";
my $lastquery = "";
my $hitnum = 0;
while(<DOT>){
   $_ =~ s/[\n\r]//g;
   my ($querygene, $sbjctgene, $identity, $matchlen, $mismatchnum, $gaplen, $querystart, $queryend, $sbjctstart, $sbjctend, $evalue, $score) = split(/\t/, $_);
   my $querychr = $querygene2chr{$querygene};
   my $sbjctchr = $sbjctgene2chr{$sbjctgene};
# only selected chromosome;
   my $is2skip1 = 1;
   for(my $i=0; $i<=$#querychr; $i++){if($querychr eq $querychr[$i]){$is2skip1 = 0; last;}}
   my $is2skip2 = 1;
   for(my $i=0; $i<=$#sbjctchr; $i++){if($sbjctchr eq $sbjctchr[$i]){$is2skip2 = 0; last;}}
   if($is2skip1 eq 1 || $is2skip2 eq 1){next;}
# max hit can not more than HITNUM;
   if($lastquery ne $querygene){$hitnum = 1;$lastquery = $querygene;}
   else{$hitnum ++;}
   if($hitnum > $_hitnum){next;}
# blast score should be more than SCORE
   if($score < $_score){next;}
# e-value should be less than EVALUE
   if($evalue > $_evalue){next;}
# get x, y value of dot
   my ($posx1, $posy1, $posx2, $posy2, $selfhit1x, $selfhit1y, $selfhit2x, $selfhit2y);
   if($poro eq "n"){
      if($querychr2order{$querychr} eq 0){$posx1 = $left_curb*1.25 + $querygene2order{$querygene}/$scale_ratio1;}
      else{$posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2order{$querygene}/$scale_ratio1;}
      if($sbjctchr2order{$sbjctchr} eq 0){$posy1 = $top_curb*1.25 + $sbjctgene2order{$sbjctgene}/$scale_ratio2;}
      else{$posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2order{$sbjctgene}/$scale_ratio2;}
    }
   elsif($poro eq "y"){
      if($querychr2order{$querychr} eq 0){$posx1 = $left_curb*1.25 + $querygene2pos{$querygene}/$scale_ratio1;}
      else{$posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2pos{$querygene}/$scale_ratio1;}
      if($sbjctchr2order{$sbjctchr} eq 0){$posy1 = $top_curb*1.25 + $sbjctgene2pos{$sbjctgene}/$scale_ratio2;}
      else{$posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2pos{$sbjctgene}/$scale_ratio2;}
    }  
# set dot's color
   my $color = $gray;
   if($hitnum == 1){$color = $red;}
   elsif($hitnum == 2){$color = $blue;}
# draw dot
   
   $img -> filledArc($posx1, $posy1, $dot_pix, $dot_pix,0,360,$color);
}
close($plot_file);
############################################################################
#########Bravo, now, finished all dotplot part drawing!
############################################################################

open(FIG,">".$output_fig) or die "can't open figure!\n";
binmode FIG;
print FIG $img -> jpeg;
close(FIG);


sub sum_chr_len(){
   my @chrlen = @_;
   my $sumlen = $chrlen[0];
   for(my $i=1; $i<=$#chrlen; $i++){
     $sumlen = $sumlen + $chrlen[$i];
   }
   return $sumlen;
}