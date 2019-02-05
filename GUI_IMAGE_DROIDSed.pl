#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";
use File::Copy;
use List::Util qw( min );
use List::Util qw( max );
use List::Util qw(min max);
use Statistics::Descriptive();


# specify the path to working directory for Chimera here
open(IN, "<"."paths.ctl") or die "could not find paths.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $path = @INrow[1];
	 if ($header eq "chimera_path"){$chimera_path = $path;}
}
close IN;
print "path to Chimera .exe\t"."$chimera_path\n";

#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- visual toolbox for functional evolutionary comparison
  of molecular dynamic simulation \n\n";

#print "Enter residue number at the start of both chains\n";
#print "(e.g. enter 389 if starts at THR 389.A) \n";
#print "(e.g. enter 1 if starts at MET 1.A) \n\n";
#my $startN = <STDIN>;
#chop($startN);

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
my $cutoffValue = 0.01;
my $colorType = '';
my $testType = '';
my $attr = '';
my $mtc = '';
my $statType = '';
my $max_val = 0;
my $min_val = 0;
my $frameCount = '';


# read control files
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      if ($header eq "reference"){$refID = $value;}
      if ($header eq "length"){$lengthID = $value;}
      if ($header eq "start"){$startN = $value;}
      if ($header eq "cutoff_value"){$cutoffValue = $value;}
      if ($header eq "test_type"){$testStr = $value;}
      if ($header eq "homology"){$homology = $value;}
      if ($header eq "num_chains"){$chainN = $value;}
}
close IN;

open(IN2, "<"."MDr.ctl") or die "could not find MDr.ctl file\n";
my @IN2 = <IN2>;
for (my $i = 0; $i < scalar @IN2; $i++){
	 my $INrow = $IN2[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "Number_Runs"){$number_runs = $value;}
      if ($header eq "Solvation_Method"){$solvation_method = $value;}
      
}
close IN2;


#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("DROIDS - Representation in UCSF Chimera"); # Titles the main window
$mw->setPalette("gray");

# Representations Frame
my $repFrame = $mw->Frame(	-label => "PROTEIN REPRESENTATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $surfaceCheck = $repFrame->Checkbutton( -text => "surface",
#						-value=>1,
                              -foreground => 'chocolate4',
						-variable=>\$surface
						);
	my $ribbonCheck = $repFrame->Checkbutton( -text => "ribbon",
#						-value=>1,
						-foreground => 'chocolate4',
                              -variable=>\$ribbon
						);


# Attribute Frame
my $attrFrame = $mw->Frame(	-label => "ATTRIBUTE TO COLOR BY",
				-relief => "groove",
				-borderwidth => 2
				);
	my $deltaRadio = $attrFrame->Radiobutton( -text => "dFLUX (as Avg difference)",
						-foreground => 'navy',
                              -value=>"delta",
						-variable=>\$attr
						);
     my $klRadio = $attrFrame->Radiobutton( -text => "dFLUX (as KL divergence)",
						-foreground => 'navy',
                              -value=>"deltaKL",
						-variable=>\$attr
						);
	my $pValueRadio = $attrFrame->Radiobutton( -text => "p-value",
						-foreground => 'navy',
                              -value=>"pval",
						-variable=>\$attr
						);
	my $statRadio = $attrFrame->Radiobutton( -text => "KS test (D stat)",
						-foreground => 'navy',
                              -value=>"dval",
						-variable=>\$attr
						);
    my $gdistRadio = $attrFrame->Radiobutton( -text => "Grantham Score",
						-foreground => 'navy',
                              -value=>"gdist",
						-variable=>\$attr
						);



# Color Type Frame
my $seqFrame = $mw->Frame(	-label => "PROTEIN COLOR SCHEME",
				-relief => "groove",
				-borderwidth => 2
				);
     my $col1Radio = $seqFrame->Radiobutton(-text=>"'stoplight' color scheme = red/green for dFLUX, blue/white for p-value, yellow for mutation",
						-foreground => 'maroon4',
                              -value=>"c1",
						-variable=>\$colorScheme
                              );
	my $col2Radio = $seqFrame->Radiobutton(-text=>"'temperature' color scheme = blue/red for dFLUX, green/white for p-value, tan for mutation",
						-foreground => 'maroon4',
                              -value=>"c2",
						-variable=>\$colorScheme
                              );
     
   
# color scale Frame
my $scaleFrame = $mw->Frame(	-label => "SCALE RANGE FOR dFLUX",
				-relief => "groove",
				-borderwidth => 2
                );
	my $autoRadio = $scaleFrame->Radiobutton( -text => "autoscale to min/max value",
						-foreground => 'darkgreen',
                              -value=>"auto",
						-variable=>\$scaleType
						);
    my $fixed2Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-2 to 2 angstrom)",
						-foreground => 'darkgreen',
                              -value=>"fixed2",
						-variable=>\$scaleType
						);
	my $fixed1Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-1 to 1 angstrom)",
						-foreground => 'darkgreen',
                              -value=>"fixed1",
						-variable=>\$scaleType
						);
    my $fixed05Radio = $scaleFrame->Radiobutton( -text => "fixed scale (-0.5 to 0.5 angstrom)",
						-foreground => 'darkgreen',
                              -value=>"fixed05",
						-variable=>\$scaleType
						);



# frameCount Frame				

my $frameFrame = $mw->Frame();
		my $frameLabel = $frameFrame->Label(-text=>"number of movie frames (e.g. 0.5ns production run = 2500 frames): ");
		my $frameEntry = $frameFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$frameCount
					);
		
# Buttons

my $displayButton = $mw -> Button(-text => "display statistics on PDB reference structure", 
				-command => \&display,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a final results display button
my $playXYZbutton = $mw -> Button(-text => "play movies on XYZ axes in DROIDS viewer", 
				-command => \&playXYZ,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a final results display button
my $playROLLbutton = $mw -> Button(-text => "play movies on ROLLING axes in DROIDS viewer", 
				-command => \&playROLL,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a final results display button
my $moviesButton = $mw -> Button(-text => "render movies in UCSF Chimera", 
				-command => \&movies,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a movie maker button

my $stopButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a exit button
my $seriesButton = $mw -> Button(-text => "create machine learning classifier on dynamics", 
				-command => \&series,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a time series button
my $attrButton = $mw -> Button(-text => "make attribute file for Chimera", 
				-command => \&attribute_file,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a attr file button
my $ctlButton = $mw -> Button(-text => "edit image control file (DROIDS.ctl)", 
				-command => \&ctl,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button

#### Organize GUI Layout ####
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$seriesButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$playXYZbutton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$playROLLbutton->pack(-side=>"bottom",
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
$ctlButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$frameFrame->pack(-side=>"bottom",
		-anchor=>"n"
		);
$frameLabel->pack(-side=>"left");
$frameEntry->pack(-side=>"left");

$surfaceCheck->pack(-anchor=>"w");
$ribbonCheck->pack(-anchor=>"w");
$repFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$col1Radio->pack(-anchor=>"w");
$col2Radio->pack(-anchor=>"w");
$deltaRadio->pack(-anchor=>"w");
$klRadio->pack(-anchor=>"w");
$pValueRadio->pack(-anchor=>"w");
$statRadio->pack(-anchor=>"w");
$gdistRadio->pack(-anchor=>"w");
$autoRadio->pack(-anchor=>"w");
$fixed2Radio->pack(-anchor=>"w");
$fixed1Radio->pack(-anchor=>"w");
$fixed05Radio->pack(-anchor=>"w");


$attrFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$seqFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$scaleFrame->pack(-side=>"top",
		-anchor=>"n"
		);
 


MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

sub stop {exit;}

########################################################################################
sub ctl {
     if ($attr eq "dval"){$colorType = "wg";}
     if ($attr eq "gdist"){$colorType = "wo";}
     if ($attr eq "delta" && $colorScheme eq "c1" ){$colorType = "rg";}
     if ($attr eq "deltaKL" && $colorScheme eq "c1" ){$colorType = "rg";}
     if ($attr eq "pval" && $colorScheme eq "c1" ){$colorType = "bw";}
     if ($attr eq "delta" && $colorScheme eq "c2" ){$colorType = "br";}
     if ($attr eq "deltaKL" && $colorScheme eq "c2" ){$colorType = "br";}
     if ($attr eq "pval" && $colorScheme eq "c2" ){$colorType = "rw";}
     if ($homology eq "loose"){$mutType = "gray50";}
     if ($homology eq "strict" && $colorScheme eq "c1"){$mutType = "yellow";}
     if ($homology eq "strict" && $colorScheme eq "c2"){$mutType = "tan";}

# make control file for DROIDS	
print("Making ctl file...\n");
	if ($ribbon == 1 && $surface == 0) {$repStr = "ribbon";}  # opaque ribbon rep only
	if ($surface == 1 && $ribbon == 0) {$repStr =  "surface";} # opaque surface rep only
	if ($surface == 1 && $ribbon == 1) {$repStr =  "ribbonsurface";} # opaque ribbon with transparent surface

	$testStr = "flux"; $testStrLong = "fluctuation";  # file and folder labels
	
open(CTL, '>', "DROIDS.ctl") or die "Could not open output file";
print CTL "query\t"."$queryID\t # Protein Data Bank ID for query structure\n";
print CTL "reference\t"."$refID\t # Protein Data Bank ID for reference structure (or neutral model)\n";
print CTL "length\t"."$lengthID\t # number of amino acids on chain\n";
print CTL "num_chains\t"."$chainN\t # number chains in PDB structure\n";
print CTL "homology\t"."$homology\t # homology as 'strict' or 'loose'\n";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "representations\t"."$repStr\t # methods of molecular representation in Chimera (ribbon and/or surface)\n";
print CTL "test_type\t"."$testStr\t # test method (sequence = local Grantham dist, structure = RMSD, fluctuation = MD)\n";
print CTL "color_scheme\t"."$colorType\t # output color scheme (stoplight=red-green, temperature=blue-red)\n";
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

if($attr eq "delta") {
	$input_file = "DROIDS$testStrLong"."AVGchain.txt";
	$relevant_column = 5;
}
if($attr eq "deltaKL") {
	$input_file = "DROIDS$testStrLong"."AVGchain.txt";
	$relevant_column = 7;
}
if($attr eq "gdist") {
	$input_file = "myGranthamDistances.txt";
	$relevant_column = 1;
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
if($attr eq "gdist") {$max_val = 215; $min_val = 0;}
if($attr eq "dval") {$max_val = $max; $min_val = 0;}
if($attr eq "pval") {$max_val = $cutoffValue; $min_val = 0;}
if($attr eq "delta" & $scaleType eq "auto") {$max_val = $max_abs_val; $min_val = -1 * $max_abs_val;}
if($attr eq "delta" & $scaleType eq "fixed2") {$max_val = 2; $min_val = -2;}
if($attr eq "delta" & $scaleType eq "fixed1") {$max_val = 1; $min_val = -1;}
if($attr eq "delta" & $scaleType eq "fixed05") {$max_val = 0.5; $min_val = -0.5;}
if($attr eq "deltaKL") {$max_val = $max_abs_val; $min_val = -1*$max_abs_val;}

sleep(1);
if (! -e "attribute_files") {mkdir "attribute_files";}
open(OUTPUT, ">"."attribute_files/$output_file");
 print OUTPUT "recipient: $unit\nattribute: $attr\n\n";
 for (my $b = 1; $b < scalar @IN; $b++) {
	print OUTPUT $OUTrows[$b]
 }
close(OUTPUT);
print("Attribute file complete\n\n");
sleep(1);
print("Renaming and copying topology and binary files for colormapping and rendering\n");
copy("wat_$refID"."REDUCED.prmtop", "wat_$refID".".prmtop");
copy("vac_$refID"."REDUCED.prmtop", "vac_$refID".".prmtop");
copy("prod_$refID"."REDUCED_0.nc", "prod_$refID"."_0.nc");
print("copying files complete\n\n");
}


############################################################################################################
sub movies{
print("Rendering 8 movies on XYZ axes...\n");
print("this may take several minutes...\n\n");
print("close Chimera window when 8 movie files appear in movies folder\n\n");
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
sub playXYZ{
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'R1', 'R2');
for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("python DROIDS_gstreamer.py @movies");
}
#############################################################################################################
sub playROLL{
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("python DROIDS_gstreamer.py @movies");
}
#############################################################################################################
sub series {

# choose topology file     
if($solvation_method eq "implicit"){$TOPfile = "vac_1ubqREDUCED.prmtop";}
if($solvation_method eq "explicit"){$TOPfile = "wat_1ubqREDUCED.prmtop";}
# loop through .nc files on refID
for(my $j = 0; $j<$number_runs; $j++){
$TRAJfile = "prod_$refID"."REDUCED_$j.nc";
$OUTfile = "fluxtime_$refID"."_$j.txt";
$step = 5;
$steplimit = $frameCount;
$start = 0;
$stop = 5;
open (CPPTRAJ, "|"."cpptraj -p $TOPfile\n");
print CPPTRAJ "trajin $TRAJfile\n";
print CPPTRAJ "rms first\n";
print CPPTRAJ "average crdset MyAvg\n";
print CPPTRAJ "run\n";
print CPPTRAJ "rms ref MyAvg\n";
print CPPTRAJ "rms first average\n";
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $steplimit\n";
print CPPTRAJ "run\n";
for(my $i = 0; $i<$steplimit/$step; $i++){
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $stop\n";
print CPPTRAJ "run\n";
$start = $start + $step;
$stop = $stop + $step;
}
print CPPTRAJ "quit\n";
close CPPTRAJ;
}
# loop through .nc files on queryID
for(my $j = 0; $j<$number_runs; $j++){
$TRAJfile = "prod_$queryID"."REDUCED_$j.nc";
$OUTfile = "fluxtime_$queryID"."_$j.txt";
$step = 5;
$steplimit = $frameCount;
$start = 0;
$stop = 5;
open (CPPTRAJ, "|"."cpptraj -p $TOPfile\n");
print CPPTRAJ "trajin $TRAJfile\n";
print CPPTRAJ "rms first\n";
print CPPTRAJ "average crdset MyAvg\n";
print CPPTRAJ "run\n";
print CPPTRAJ "rms ref MyAvg\n";
print CPPTRAJ "rms first average\n";
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $steplimit\n";
print CPPTRAJ "run\n";
for(my $i = 0; $i<$steplimit/$step; $i++){
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $stop\n";
print CPPTRAJ "run\n";
$start = $start + $step;
$stop = $stop + $step;
}
print CPPTRAJ "quit\n";
close CPPTRAJ;
}
sleep(1); print("time series files created\n\n"); sleep(1);    
     
}
###########################################################################################################