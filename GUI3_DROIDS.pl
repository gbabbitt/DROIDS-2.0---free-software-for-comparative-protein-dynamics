#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";
use File::Copy;
use List::Util qw( min );
use List::Util qw( max );
use List::Util qw(min max);
use Descriptive();

# specify the path to working directory for Chimera here

my $chimera_path = "/opt/UCSF/Chimera64-1.11/bin/";


#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- visual toolbox for functional evolutionary comparison
  of molecular dynamic simulation \n\n";

print "Enter residue number at the start of both chains\n";
print "(e.g. enter 389 if starts at THR 389.A) \n";
print "(e.g. enter 1 if starts at MET 1.A) \n\n";
my $startN = <STDIN>;
chop($startN);

#### Declare variables ####
my @rep = ();
my @test = ();
my $queryID = '';
my $refID = '';
my $lengthID = '';
my $repStr = '';
my $testStr = '';
my $surface = 0;
my $ribbon = 0;
my $rep;
#my $flux = 0;
#my $corr = 0;
my $cutoffValue = 0.05;
my $colorType = '';
my $testType = '';
my $attr = '';
my $mtc = '';
my $statType = '';
my $max_val = 0;
my $min_val = 0;
my $frameCount = '';

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("Welcome to DROIDS!"); # Titles the main window


my $cutoffScale = $mw->Scale(-label=>"Cut-off for Significance :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>0.1,
			-variable=>\$cutoffValue,
			-tickinterval=>0.05,
			-resolution=>0.001,
			-length=>200,
			);


# Representations Frame
my $repFrame = $mw->Frame(	-label => "PROTEIN REPRESENTATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $surfaceCheck = $repFrame->Checkbutton( -text => "surface",
#						-value=>1,
						-variable=>\$surface
						);
	my $ribbonCheck = $repFrame->Checkbutton( -text => "ribbon",
#						-value=>1,
						-variable=>\$ribbon
						);


# Color Type Frame
my $seqFrame = $mw->Frame(	-label => "PROTEIN COLOR SCHEME",
				-relief => "groove",
				-borderwidth => 2
				);
	my $rgRadio = $seqFrame->Radiobutton(-text=>"red-green (for delta)",
						-value=>"rg",
						-variable=>\$colorType
						);
	my $ybRadio = $seqFrame->Radiobutton(-text=>"yellow-blue (for delta)",
						-value=>"yb",
						-variable=>\$colorType
						);
	my $omRadio = $seqFrame->Radiobutton(-text=>"orange-magenta (for delta)",
						-value=>"om",
						-variable=>\$colorType
						);	
	my $bwRadio = $seqFrame->Radiobutton(-text=>"blue-white (for p-value)",
						-value=>"bw",
						-variable=>\$colorType
						);
	my $rwRadio = $seqFrame->Radiobutton(-text=>"red-white (for p-value)",
						-value=>"rw",
						-variable=>\$colorType
						);
	my $wgRadio = $seqFrame->Radiobutton(-text=>"white-green (for test stat)",
						-value=>"wg",
						-variable=>\$colorType
						);
	my $woRadio = $seqFrame->Radiobutton(-text=>"white-orange (for test stat)",
						-value=>"wo",
						-variable=>\$colorType
						);
	
# Test Type Frame
my $testFrame = $mw->Frame(	-label => "ATOM MOTION TYPE TO ANALYZE",
				-relief => "groove",
				-borderwidth => 2
				);
	my $fluxCheck = $testFrame->Radiobutton( -text => "atom fluctuation",
						-value=>"fx",
						-variable=>\$testType
						);
	my $corrCheck = $testFrame->Radiobutton( -text => "atom correlation",
						-value=>"cr",
						-variable=>\$testType
						);

# Attribute Frame
my $attrFrame = $mw->Frame(	-label => "ATTRIBUTE TO COLOR BY",
				-relief => "groove",
				-borderwidth => 2
				);
	my $deltaRadio = $attrFrame->Radiobutton( -text => "delta (flux/corr)",
						-value=>"delta",
						-variable=>\$attr
						);
	my $pValueRadio = $attrFrame->Radiobutton( -text => "p-value",
						-value=>"pval",
						-variable=>\$attr
						);
	my $statRadio = $attrFrame->Radiobutton( -text => "test stat (D value)",
						-value=>"dval",
						-variable=>\$attr
						);

# Scaling Frame
my $scalingFrame = $mw->Frame(	-label => "dFLUX CALCULATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $absoluteRadio = $scalingFrame->Radiobutton( -text => "absolute   (use unmodified FLUX values)",
						-value=>"absolute",
						-variable=>\$scalingType
						);
	my $relativeRadio = $scalingFrame->Radiobutton( -text => "relative   (scale query FLUX to avg reference FLUX)",
						-value=>"relative",
						-variable=>\$scalingType
						);


# multiple test correction Frame
my $mtcFrame = $mw->Frame(	-label => "MULTIPLE TEST CORRECTION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $noneRadio = $mtcFrame->Radiobutton( -text => "none   (not recommended)",
						-value=>"none",
						-variable=>\$mtc
						);
	my $bonfRadio = $mtcFrame->Radiobutton( -text => "Bonferroni   (can be overly strict)",
						-value=>"bonferroni",
						-variable=>\$mtc
						);
	my $fdrRadio = $mtcFrame->Radiobutton( -text => "Benjamini-Hochberg   (recommended)",
						-value=>"fdr",
						-variable=>\$mtc
						);

# color scale Frame
my $scaleFrame = $mw->Frame(	-label => "SCALE RANGE for\n delta(flux/corr)\n",
				-relief => "groove",
				-borderwidth => 2
				);
	my $autoRadio = $scaleFrame->Radiobutton( -text => "autoscale to min/max value",
						-value=>"auto",
						-variable=>\$scaleType
						);
    my $fixed2Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-2 to 2 angstrom)",
						-value=>"fixed2",
						-variable=>\$scaleType
						);
	my $fixed1Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-1 to 1 angstrom)",
						-value=>"fixed1",
						-variable=>\$scaleType
						);
    my $fixed05Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-0.5 to 0.5 angstrom)",
						-value=>"fixed05",
						-variable=>\$scaleType
						);


# mutation color Frame
my $mutFrame = $mw->Frame(	-label => "COLOR - NONHOMOLOGOUS REGIONS & MUTATIONS)",
				-relief => "groove",
				-borderwidth => 2
				);
	my $yellowRadio = $mutFrame->Radiobutton( -text => "yellow - (strict homology - small distance)",
						-value=>"yellow",
						-variable=>\$mutType
						);
    my $redRadio = $mutFrame->Radiobutton( -text => "red -(strict homology - small distance)",
						-value=>"red",
						-variable=>\$mutType
						);
	my $tanRadio = $mutFrame->Radiobutton( -text => "tan - (loose homology - large distance)",
						-value=>"tan",
						-variable=>\$mutType
						);
    my $grayRadio = $mutFrame->Radiobutton( -text => "dark gray - (loose homology - large distance)",
						-value=>"gray50",
						-variable=>\$mutType
						);


# PDB ID Frame				
my $pdbFrame = $mw->Frame();
	my $queryFrame = $pdbFrame->Frame();
		my $queryLabel = $queryFrame->Label(-text=>"pdb ID query (e.g. 1ubq): ");
		my $queryEntry = $queryFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$queryID
					);
	my $refFrame = $pdbFrame->Frame();
		my $refLabel = $refFrame->Label(-text=>"pdb ID reference (e.g. 2ubq): ");
		my $refEntry = $refFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$refID
					);
	my $lengthFrame = $pdbFrame->Frame();
		my $lengthLabel = $lengthFrame->Label(-text=>"length of protein (no. of AA's): ");
		my $lengthEntry = $lengthFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$lengthID
					);
	my $frameFrame = $pdbFrame->Frame();
		my $frameLabel = $frameFrame->Label(-text=>"length of movies (e.g. 500 frames): ");
		my $frameEntry = $frameFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$frameCount
					);	
		
		
		
# Buttons

my $statsButton = $mw -> Button(-text => "make statistical comparisons and plots in R", 
				-command => \&stats
				); # Creates a stats test button

my $displayButton = $mw -> Button(-text => "display statistics on PDB reference structure", 
				-command => \&display
				); # Creates a final results display button
my $playButton = $mw -> Button(-text => "play movies on XYZ axes in DROIDS viewer", 
				-command => \&play
				); # Creates a final results display button

my $moviesButton = $mw -> Button(-text => "render movies on XYZ axes in Chimera", 
				-command => \&movies
				); # Creates a movie maker button

my $stopButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop
				); # Creates a exit button

my $ctlButton = $mw -> Button(-text => "make control file (DROIDS.ctl)", 
				-command => \&ctl
				); # Creates a ctl file button

my $attrButton = $mw -> Button(-text => "make attribute file for Chimera", 
				-command => \&attribute_file
				); # Creates a attr file button


#### Organize GUI Layout ####
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$playButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$moviesButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$displayButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$attrButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$statsButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$ctlButton->pack(-side=>"bottom",
			-anchor=>"s"
			);


$queryLabel->pack(-side=>"left");
$queryEntry->pack(-side=>"left");
$refLabel->pack(-side=>"left");
$refEntry->pack(-side=>"left");
$lengthLabel->pack(-side=>"left");
$lengthEntry->pack(-side=>"left");
$frameLabel->pack(-side=>"left");
$frameEntry->pack(-side=>"left");
$queryFrame->pack(-side=>"top",
		-anchor=>"e");
$refFrame->pack(-side=>"top",
		-anchor=>"e");
$lengthFrame->pack(-side=>"top",
		-anchor=>"e");
$frameFrame->pack(-side=>"top",
		-anchor=>"e");
$pdbFrame->pack(-side=>"top",
		-anchor=>"n");
$testFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$cutoffScale->pack(-side=>"top");
$mtcFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$scalingFrame->pack(-side=>"top",
		-anchor=>"s"
		);

$surfaceCheck->pack();
$ribbonCheck->pack();
$fluxCheck->pack();
$corrCheck->pack();
$repFrame->pack(-side=>"left",
		-anchor=>"n"
		);
$rgRadio->pack();
$ybRadio->pack();
$omRadio->pack();
$bwRadio->pack();
$rwRadio->pack();
$wgRadio->pack();
$woRadio->pack();
$deltaRadio->pack();
$pValueRadio->pack();
$statRadio->pack();
$noneRadio->pack();
$bonfRadio->pack();
$fdrRadio->pack();
$autoRadio->pack();
$fixed2Radio->pack();
$fixed1Radio->pack();
$fixed05Radio->pack();
$absoluteRadio->pack();
$relativeRadio->pack();
$yellowRadio->pack();
$redRadio->pack();
$tanRadio->pack();
$grayRadio->pack();


$seqFrame->pack(-side=>"left",
		-anchor=>"n"
		);
$attrFrame->pack(-side=>"left",
		-anchor=>"n"
		);
$scaleFrame->pack(-side=>"left",
		-anchor=>"n"
		);
$mutFrame->pack(-side=>"left",
		-anchor=>"n"
		);

MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

sub stop {exit;}

########################################################################################
sub ctl {
# make control file for DROIDS	
print("Making ctl file...\n");
	if ($ribbon == 1 && $surface == 0) {$repStr = "ribbon";}  # opaque ribbon rep only
	if ($surface == 1 && $ribbon == 0) {$repStr =  "surface";} # opaque surface rep only
	if ($surface == 1 && $ribbon == 1) {$repStr =  "ribbonsurface";} # opaque ribbon with transparent surface

	if ($testType  eq "fx") {$testStr = "flux"; $testStrLong = "fluctuation";}
	if ($testType eq "cr") {$testStr = "corr"; $testStrLong = "correlation";}
		
open(CTL, '>', "DROIDS.ctl") or die "Could not open output file";
print CTL "query\t"."$queryID\t # Protein Data Bank ID for query structure\n";
print CTL "reference\t"."$refID\t # Protein Data Bank ID for reference structure (or neutral model)\n";
print CTL "length\t"."$lengthID\t # number of amino acids on chain\n";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "representations\t"."$repStr\t # methods of molecular representation in Chimera (ribbon and/or surface)\n";
print CTL "test_type\t"."$testStr\t # methods of molecular representation in Chimera (ribbon and/or surface)\n";
print CTL "color_scheme\t"."$colorType\t # output color scheme (red-green, yellow-blue, or orange-magenta)\n";
close CTL;
print("CTL file made\n");
}

########################################################################
sub attribute_file {

# Make Chimera-readable attribute file needed
my $relevant_column = 0;
$input_folder = "DROIDS_results_$queryID"."_$refID"."_$testStr"."_$cutoffValue";

if($attr eq "pval") {
	$input_file = "adjustKStests.txt";
	$relevant_column = 4;
}
if($attr eq "dval") {
	$input_file = "adjustKStests.txt";
	$relevant_column = 3;
}

if($attr eq "delta" && $scalingType eq "relative") {
	$input_file = "DROIDS$testStrLong"."AVGscaled.txt";
	$relevant_column = 5;
}
if($attr eq "delta" && $scalingType eq "absolute") {
	$input_file = "DROIDS$testStrLong"."AVG.txt";
	$relevant_column = 5;
}
#print("$relevant_column\n");
#print "$input_folder\n";
#print "$input_file\n";
$output_file = "attr_$queryID"."_$refID"."_$attr"."_$testStr"."_$cutoffValue.dat";
$unit = "residues";
my $symbol;
if ($unit eq "residues") { $symbol = "\:"; } 
if ($unit eq "atoms") { $symbol = "\@"; }
my @valuesList = "";
open(INPUT, "<"."$input_folder/$input_file") or die "Could not find $input_file";
 my @IN = <INPUT>;
 my @OUTrows;
 for (my $a = 1; $a < scalar @IN; $a++) {
	
    my $INrow = $IN[$a];
	my @INrow = split (/\s+/, $INrow);
	my $index = $INrow[0] - ($startN - 1);
	my $value = $INrow[$relevant_column];
	#print $value;
	$valuesList[$a] = $value;
	$index = int $index;
#	if ($attr eq "pval" and $value > $cutoffValue) {
#		$OUTrows[$a] = "\t$symbol"."$index\tns\n";
#	} else { 
		$OUTrows[$a] = "\t$symbol"."$index\t$value\n";
#	}
 }
close(INPUT);
#print("@valuesList\n");
my $min = min(@valuesList);
my $max = max(@valuesList);
#$statSCORE = new Statistics::Descriptive::Full; # calc min and max value
#$statSCORE->add_data (@valuesList);
#$min = $statSCORE->min();
#$max = $statSCORE->max();
#print("max is $max\n");
#print("min is $min\n");
my $max_abs_val = abs(max($max,-$min));
#print("abs max is $max_abs_val\n");
if($attr eq "dval") {$max_val = $max; $min_val = 0;}
if($attr eq "pval") {$max_val = $cutoffValue; $min_val = 0;}
if($attr eq "delta" & $scaleType eq "auto") {$max_val = $max_abs_val; $min_val = -1 * $max_abs_val;}
if($attr eq "delta" & $scaleType eq "fixed2") {$max_val = 2; $min_val = -2;}
if($attr eq "delta" & $scaleType eq "fixed1") {$max_val = 1; $min_val = -1;}
if($attr eq "delta" & $scaleType eq "fixed05") {$max_val = 0.5; $min_val = -0.5;}

sleep(1);
if (! -e "attribute_files") {mkdir "attribute_files";}
open(OUTPUT, ">"."attribute_files/$output_file");
 print OUTPUT "recipient: $unit\nattribute: $attr\n\n";
 for (my $b = 1; $b < scalar @IN; $b++) {
	print OUTPUT $OUTrows[$b]
 }
close(OUTPUT);
print("Attribute file complete\n");
}


#################################################################

sub stats {

# run KS tests
print " running KS tests on each amino acid\n\n";
sleep(1);
print " reading control file\n\n";

my $queryID = '';
my $referenceID = '';
my $AA_count = '';
my $level_sig = '';
my $surface_or_ribbon = '';
my $flux_or_corr = '';
my $color_scheme = '';


open(IN, "<"."DROIDS.ctl") or die "could not find CPPTRAJ input control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "query") { $queryID = $value;}
    if ($header eq "reference") { $referenceID = $value;}
    if ($header eq "length") { $AA_count = $value;}
    if ($header eq "cutoff_value") { $level_sig = $value;}
    if ($header eq "representations") { $surface_or_ribbon = $value;}
    if ($header eq "test_type") { $flux_or_corr = $value;}
	if ($header eq "color_scheme") { $color_scheme = $value;}

}
close IN;
sleep(1);
##########################################
if ($flux_or_corr eq "flux" && $scalingType eq "relative"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
    
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomfluxscaled/DROIDSfluctuation_$filenumber.txt") or next;
    my @INFO = <INFO>;
    my $INFOrow = $INFO[1];
	             my @INFOrow = split(/\s+/, $INFOrow); 
			     my $posAA = $INFOrow[1];
				 my $refAA = $INFOrow[2];
				 my $queryAA = $INFOrow[3];
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   my $sample = ''; # sample number
   my $pos_ref = ''; # amino acid position
   my $res_ref = ''; # reference residue
   my $res_query = ''; # query residue
   my $atomnumber = ''; # cpptraj atom number
   my $atomlabel = ''; # cpptraj atom number
   my $flux_ref = ''; # flux on reference residue
   my $flux_query = ''; # flux on query residue
    # read data into R
   print Rinput "data = read.table('atomfluxscaled/DROIDSfluctuation_$filenumber.txt', header = TRUE)\n"; 
   $sample = "data\$sample"; # sample number
   $pos_ref = "data\$pos_ref"; # amino acid position
   $res_ref = "data\$res_ref"; # reference residue
   $res_query = "data\$res_query"; # query residue
   $atomnumber = "data\$atom number"; # cpptraj atom number
   $atomlabel = "data\$atom label"; # cpptraj atom number
   $flux_ref = "data\$flux_ref"; # flux on reference residue
   $flux_query = "data\$flux_query"; # flux on query residue
   print Rinput "d1 = data.frame(pos=$pos_ref, res=$res_ref, fluxR=$flux_ref, fluxQ=$flux_query)\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   #print to file
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
   my @IN2 = <IN2>;
   my $label_sig = '';
   for (my $rr = 0; $rr < scalar @IN2; $rr++){ # scan atom type
			     my $IN2row = $IN2[$rr];
	             my @IN2row = split(/\s+/, $IN2row); 
			     my $testD = $IN2row[0];
				 my $Dval = $IN2row[2];
				 $Dval =~ s/,//g; # remove trailing comma
				 $Dval = $Dval + 0; # force to a number 
				 my $pval = $IN2row[5];
				 $pval = $pval + 0; # force to a number 
				 if ($pval <= $level_sig){$label_sig = "<1alpha";}
				 if ($pval <= (0.5*$level_sig)){$label_sig = "<0.5alpha";}
				 if ($pval <= (0.1*$level_sig)){$label_sig = "<0.1alpha";}
				 if ($pval > $level_sig){$label_sig = "ns";}
				 if ($testD eq "D"){print STAT2 "$posAA\t"."$refAA\t"."$queryAA\t"."$Dval\t"."$pval\t"."$label_sig\n"}
                 }
   close IN2;
   
   }	
close STAT2;
## remove temp file
unlink("./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVG.txt", "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVG.txt");
copy("DROIDSfluctuationAVGscaled.txt", "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVGscaled.txt");

}

##########################################
if ($flux_or_corr eq "flux" && $scalingType eq "absolute"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomflux/DROIDSfluctuation_$filenumber.txt") or next;
    my @INFO = <INFO>;
    my $INFOrow = $INFO[1];
	             my @INFOrow = split(/\s+/, $INFOrow); 
			     my $posAA = $INFOrow[1];
				 my $refAA = $INFOrow[2];
				 my $queryAA = $INFOrow[3];
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   my $sample = ''; # sample number
   my $pos_ref = ''; # amino acid position
   my $res_ref = ''; # reference residue
   my $res_query = ''; # query residue
   my $atomnumber = ''; # cpptraj atom number
   my $atomlabel = ''; # cpptraj atom number
   my $flux_ref = ''; # flux on reference residue
   my $flux_query = ''; # flux on query residue
    # read data into R
   print Rinput "data = read.table('atomflux/DROIDSfluctuation_$filenumber.txt', header = TRUE)\n"; 
   $sample = "data\$sample"; # sample number
   $pos_ref = "data\$pos_ref"; # amino acid position
   $res_ref = "data\$res_ref"; # reference residue
   $res_query = "data\$res_query"; # query residue
   $atomnumber = "data\$atom number"; # cpptraj atom number
   $atomlabel = "data\$atom label"; # cpptraj atom number
   $flux_ref = "data\$flux_ref"; # flux on reference residue
   $flux_query = "data\$flux_query"; # flux on query residue
   print Rinput "d1 = data.frame(pos=$pos_ref, res=$res_ref, fluxR=$flux_ref, fluxQ=$flux_query)\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   #print to file
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
   my @IN2 = <IN2>;
   my $label_sig = '';
   for (my $rr = 0; $rr < scalar @IN2; $rr++){ # scan atom type
			     my $IN2row = $IN2[$rr];
	             my @IN2row = split(/\s+/, $IN2row); 
			     my $testD = $IN2row[0];
				 my $Dval = $IN2row[2];
				 $Dval =~ s/,//g; # remove trailing comma
				 $Dval = $Dval + 0; # force to a number 
				 my $pval = $IN2row[5];
				 $pval = $pval + 0; # force to a number 
				 if ($pval <= $level_sig){$label_sig = "<1alpha";}
				 if ($pval <= (0.5*$level_sig)){$label_sig = "<0.5alpha";}
				 if ($pval <= (0.1*$level_sig)){$label_sig = "<0.1alpha";}
				 if ($pval > $level_sig){$label_sig = "ns";}
				 if ($testD eq "D"){print STAT2 "$posAA\t"."$refAA\t"."$queryAA\t"."$Dval\t"."$pval\t"."$label_sig\n"}
                 }
   close IN2;
   
   }	
close STAT2;
## remove temp file
unlink("./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVG.txt", "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVG.txt");
copy("DROIDSfluctuationAVGscaled.txt", "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVGscaled.txt");

}

########################################
if ($flux_or_corr eq "corr"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomcorr/DROIDScorrelation_$filenumber.txt") or next;
    my @INFO = <INFO>;
    my $INFOrow = $INFO[1];
	             my @INFOrow = split(/\s+/, $INFOrow); 
			     my $posAA = $INFOrow[1];
				 my $refAA = $INFOrow[2];
				 my $queryAA = $INFOrow[3];
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   my $sample = ''; # sample number
   my $pos_ref = ''; # amino acid position
   my $res_ref = ''; # reference residue
   my $res_query = ''; # query residue
   my $atomnumber = ''; # cpptraj atom number
   my $atomlabel = ''; # cpptraj atom number
   my $corr_ref = ''; # corr on reference residue
   my $corr_query = ''; # corr on query residue

    # read data into R
   print Rinput "data = read.table('atomcorr/DROIDScorrelation_$r.txt', header = TRUE)\n"; 
   $sample = "data\$sample"; # sample number
   $pos_ref = "data\$pos_ref"; # amino acid position
   $res_ref = "data\$res_ref"; # reference residue
   $res_query = "data\$res_query"; # query residue
   $atomnumber = "data\$atom number"; # cpptraj atom number
   $atomlabel = "data\$atom label"; # cpptraj atom number
   $corr_ref = "data\$corr_ref"; # flux on reference residue
   $corr_query = "data\$corr_query"; # flux on query residue
   
   print Rinput "d1 = data.frame(pos=$pos_ref, res=$res_ref, corrR=$corr_ref, corrQ=$corr_query)\n";
   print Rinput "ks_test<-ks.test($corr_ref, $corr_query)\n";
   print Rinput "print(ks_test)\n";
   #print to file
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($corr_ref, $corr_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
   my @IN2 = <IN2>;
   my $label_sig = '';
   for (my $rr = 0; $rr < scalar @IN2; $rr++){ # scan atom type
			     my $IN2row = $IN2[$rr];
	             my @IN2row = split(/\s+/, $IN2row); 
			     my $testD = $IN2row[0];
				 my $Dval = $IN2row[2];
				 $Dval =~ s/,//g; # remove trailing comma
				 $Dval = $Dval + 0; # force to a number 
				 my $pval = $IN2row[5];
				 $pval = $pval + 0; # force to a number 
				 if ($pval <= $level_sig){$label_sig = "<1alpha";}
				 if ($pval <= (0.5*$level_sig)){$label_sig = "<0.5alpha";}
				 if ($pval <= (0.1*$level_sig)){$label_sig = "<0.1alpha";}
				 if ($pval > $level_sig){$label_sig = "ns";}
				 if ($testD eq "D"){print STAT2 "$posAA\t"."$refAA\t"."$queryAA\t"."$Dval\t"."$pval\t"."$label_sig\n"}
                 }
   close IN2;
   
   }	
close STAT2;
## remove temp file
unlink("./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDScorrelationAVG.txt", "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDScorrelationAVG.txt");
}

################## adjust KStests.txt for multiple tests ############

sleep(1);
print " calculating multiple test correction \n\n";
sleep(1);
open (OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/adjustKStests.txt") or die "could not create output file\n";
print OUT "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";
open (TMP, ">"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStests.txt', header = TRUE)\n"; 
$posAA = "data\$posAA"; # position on reference structure	
$refAA = "data\$refAA"; # AA label on reference structure
$queryAA = "data\$queryAA"; # AA label on query structure
$Dval = "data\$Dval"; # D value for KS test
$pval = "data\$pval"; # p value for KS test
$signif = "data\$signif"; # significance label
print Rinput "p.adjust($pval, method = 'bonferroni', n = length($pval))\n";
print Rinput "p.adjust($pval, method = 'fdr', n = length($pval))\n"; #adjust p values for false discovery rate (i.e. Benjamini-Hochberg procedure)
print Rinput "write(p.adjust($pval, method = '$mtc', n = length($pval)), './DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/adjustPvalues.txt', sep='\n')\n";
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";	
close TMP;
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/KStests.txt") or die "could not create output file\n";
open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
my @IN1 = <IN1>;
my @IN2 = <IN2>;
my $label_sig = '';
                 
    for (my $k = 0; $k < scalar @IN1; $k++){ # scan KStest file
			     my $IN1row = $IN1[$k];
	             my @IN1row = split(/\s+/, $IN1row); 
			     my $posAA = $IN1row[0];
				 my $refAA = $IN1row[1];
				 my $queryAA = $IN1row[2];
				 my $Dval = $IN1row[3];
				 my $oldPval = $IN1row[4];
				 my $newPval = '';
				 for (my $kk = 0; $kk < scalar @IN2; $kk++){ # scan adjusted p values
				   my $IN2row = $IN2[$kk]; my @IN2row = split(/\s+/, $IN2row);
				   my $findvalue = $IN2row[0];
				   if ($k == $kk+1){$newPval = $findvalue;}
				   }
				   if ($newPval <= $level_sig){$label_sig = "<1alpha";}
				   if ($newPval <= (0.5*$level_sig)){$label_sig = "<0.5alpha";}
				   if ($newPval <= (0.1*$level_sig)){$label_sig = "<0.1alpha";}
				   if ($newPval > $level_sig){$label_sig = "ns";}
                   if ($posAA >= 1){print OUT "$posAA\t"."$refAA\t"."$queryAA\t"."$Dval\t"."$newPval\t"."$label_sig\n";}
				   
	
}
close OUT;

################## plotting R graphics ##############################

sleep(1);
print " plotting KS tests\n\n";
sleep(1);

open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";

# read data into R
if ($flux_or_corr eq "flux" && $scalingType eq "relative"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVGscaled.txt', header = TRUE)\n";} 
if ($flux_or_corr eq "flux" && $scalingType eq "absolute"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSfluctuationAVG.txt', header = TRUE)\n";} 
if ($flux_or_corr eq "corr"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDScorrelationAVG.txt', header = TRUE)\n";} 
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/adjustKStests.txt', header = TRUE)\n"; 
if ($flux_or_corr eq "flux"){
$pos_ref = "data1\$pos_ref"; # position on reference structure	
$res_ref = "data1\$res_ref"; # AA label on reference structure
$res_query = "data1\$res_query"; # AA label on query structure
$flux_ref_avg = "data1\$flux_ref_avg"; # avg flux on reference structure
$flux_query_avg = "data1\$flux_query_avg"; # avg flux on query structure
$delta_flux = "data1\$delta_flux"; # signed difference in avg flux
$abs_delta_flux = "data1\$abs_delta_flux"; # unsigned difference in avg flux
}
if ($flux_or_corr eq "corr"){
$pos_ref = "data1\$pos_ref"; # position on reference structure	
$res_ref = "data1\$res_ref"; # AA label on reference structure
$res_query = "data1\$res_query"; # AA label on query structure
$corr_ref_avg = "data1\$corr_ref_avg"; # avg flux on reference structure
$corr_query_avg = "data1\$corr_query_avg"; # avg flux on query structure
$delta_corr = "data1\$delta_corr"; # signed difference in avg flux
$abs_delta_corr = "data1\$abs_delta_corr"; # unsigned difference in avg flux
}
$posAA = "data2\$posAA"; # position on reference structure	
$refAA = "data2\$refAA"; # AA label on reference structure
$queryAA = "data2\$queryAA"; # AA label on query structure
$Dval = "data2\$Dval"; # D value for KS test
$pval = "data2\$pval"; # p value for KS test
$signif = "data2\$signif"; # significance label
print Rinput "data1\n";
print Rinput "data2\n";
# barplot
if ($flux_or_corr eq "flux"){
print Rinput "d1A = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = 'query PDB')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = 'ref PDB')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'position (residue number)', y = 'direction dFLUX(signed)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dFLUX(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'))\n";
}

if ($flux_or_corr eq "corr"){
print Rinput "d1A = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$corr_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$corr_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref, label2r=$res_ref, label2q=$res_query, Y2val=$delta_corr)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_corr)\n";
print Rinput "d4 = data.frame(pos4=$posAA, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number) ref=red query=blue', y = 'avg CORR (ref vs query)') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = 'query PDB')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = 'ref PDB')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'position (residue number)', y = 'direction dCORR(signed)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dCORR(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'))\n";
}

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
#print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/"."DROIDSplot.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVG.txt") or die "could not open file\n";
my @IN = <IN>;
for (my $c = 0; $c < scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $posit = $INrow[0];
    my $val = $INrow[5];
    if ($posit ne "pos_ref"){push(@vals, $val);}
    $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
    $statSCORE->add_data (@vals);
	$scaling_factor = $statSCORE->mean();
}
close IN;
print "scaling factor used on query FLUX = "."$scaling_factor"." angstroms\n\n";
sleep(1);
}
if ($scalingType eq "absolute"){print "no scaling factor used\n\n"; sleep(1);}
print " stats and plotting subroutine is complete\n\n";
print " close PDF viewer to continue\n\n";
system "evince ./DROIDS_results_$queryID"."_$refID"."_$flux_or_corr"."_$level_sig"."/DROIDSplot.pdf\n";
}

############################################################################################################
sub movies{
print("Rendering 6 movies on XYZ axes...\n");
print("this may take several minutes...\n\n");
print("close Chimera window when 6 movie files appear in movies folder\n\n");
#print("/opt/UCSF/Chimera64-1.11/bin/chimera --script \"render_movies.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"render_movies_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"render_movies_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "red"){system("$chimera_path"."chimera --script \"render_movies_red.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "yellow"){system("$chimera_path"."chimera --script \"render_movies_yellow.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Movies rendered\n");
}

############################################################################################################

sub display{
print("Preparing static display...\n");
print("close Chimera window to exit\n\n");
print("ignore MD movie window unless you want to make a custom movie\n\n");
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "red"){system("$chimera_path"."chimera --script \"color_by_attr_red.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "yellow"){system("$chimera_path"."chimera --script \"color_by_attr_yellow.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");

}

############################################################################################################
sub play{
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'R1', 'R2');
for (my $i = 0; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("python DROIDS_gstreamer.py @movies");
}
#############################################################################################################