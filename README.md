# GDotplotter
A drawing dotplot software for showing genome homology
---------------------------------------------------------------------------------------------------------------------------------

Download

Source code, binaries, Perl scripts, example data, and documentation for GDotplotter can be freely obtained from the website https://github.com/tianyuwuai/GDotplotter.

Source code: GDotplotter.Source.zip

Windows release: GDotplotter.Release.zip

Perl scripts for Linux, MacOS and Windows: GDotplotter.Perl.Script.zip 

Prerequisites

For running GDotplotter’s Perl script successfully on any systems, it is needed to install the Perl and the GD module into your system. 
For running GDotplotter’s graphical user interface on Windows, ActivePerl is required and also need add Perl to PATH environment variable.

Usage

1.	blast2dotplot.pl
$ perl blast2dotplot.pl querychr sbjctchr querygff1 objctgff2 blastm8 pos_or_order queryinitial sbjctinitial e-value score hitnumber    

2.	blast2dotplot-cdsvschr.pl
$ perl blast2dotplot-cdsvschr.pl querychr sbjctchr querygff1 objctgff2 blastm8 queryinitial sbjctinitial e-value score hitnumber                 

3.	mcscan2dotplot.pl
$ perl mcscan2dotplot.pl querychr sbjctchr querygff1 objctgff2 mcscan_align pos_or_order queryinitial sbjctinitial score e-value length             

4.	adhore2dotplot.pl
$ perl adhore2dotplot.pl querychr sbjctchr querygff1 objctgff2 inputfile pos_or_order queryinitial sbjctinitial                                                                       
5.	colinearscan2dotplot.pl
$ perl colinearscan2dotplot.pl querychr sbjctchr querygff1 objctgff2 blockfile pos_or_order queryinitial sbjctinitial p-value length                     

6.	colinearscan2dotplot-genevschr.pl
$ perl colinearscan2dotplot-genevschr.pl querychr sbjctchr querygff1 objctgff2 blockfile queryinitial sbjctinitial p-value length                  

7.	comparedotplot.pl
$ perl comparedotplot.pl querychr sbjctchr querygff1 objctgff2 pos_or_order queryinitial sbjctinitial file1 ? file2 ?                                                            
Inputs Files and formats

1.	GDotplotter allows users to provide 6 type files as the toolkit input file.
a.	Blast M8 output file which nucleotide(protein) blast nucleotide(protein);
b.	Blast M8 output file which CDS blast chromosome genome;
Sb10G01042	chr06	83.33	84	14	0	286	369	7926555	7926472	2.00E-06	56
c.	MCScan .aligns file corresponding to pairwise synteny;
d.	i-ADHoRe’s anchor points table which is written into a file called anchorpoints.txt;
e.	ColinearScan block file that CDS blast chromosome sequence result M8 file as the ColinearScan input file;
f.	ColinearScan block file that nucleotide(protein) blast nucleotide(protein) result M8 file as the ColinearScan input file.

2.  GDotplotter also allows users to provide any two type files of a, c, d and f for comparing different collinear block output files which are got from different programs. 

3.  In addition, both species’ .gff files are also needed to offer each gene physical position or gene order on a chromosome in drawing a dotplot map. The .gff file contains the following tab-delimited format.
Chr01	Sb01G00001	1619	2809
Parameters
Parameterized GDotplotter is convenient for drawing specific dotplots. Peculiarly, there are some different parameters after users selecting one input file type. 
querychr	　	Query species chromosome number. Join the specific chromosome number in "_"
sbjctchr	　	Subject species chromosome number. Join the specific chromosome number in "_"
querygff1	　	Query species .gff file.
objctgff2	　	Subject species .gff file.
inputfile	　	The path of input file.
pos_or_order	　	[pos|order] Drawing by the physical location of a gene or gene order along chromosomes.
queryinitial	　	Initial of query species for headline of x axis. 
sbjctinitial	　	Initial of subject species for headline of y axis. 
Blast	　	　
e-value	　	Filter BLAST file by e-value, save all gene pairs which e-value less than setting.
score	　	Filter BLAST file by score, save all gene pairs which score more than setting.
hitnumber	　	Filter BLAST file by utmost matching gene pairs.
MCScan	　	　
score	　	Filter MCScan file by score, save blocks which score more than setting.
e-value	　	Filter MCScan file by e-value, save blocks which e-value less than setting.
length	　	Saving blocks which length more than setting.
ColinearScan	　	　
p-value	　	Filter ColinearScan file by p-value, saving smaller.
length	　	Saving blocks which length more than setting.
Example
There is a shell named “run.example.sh” to draw different dotplots by using some example data. You can just simply run it or check the details.
$ sh run.example.sh                                                                  
Changelog
Mar 26, 2016 initial release.
Contact
Any question, problem, bugs are welcome and should be dumped to
Tianyu Lei: leitianyu.china@gmail.com
