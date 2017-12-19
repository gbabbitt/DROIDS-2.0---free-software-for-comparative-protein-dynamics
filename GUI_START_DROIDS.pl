#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";
# specify the path to working directory for Chimera here
my $chimera_path = "/opt/UCSF/Chimera64-1.11/bin/";

#### Introductory Message #############

print "\n\nWelcome to DROIDS - a pipeline for evolutionary and functional comparison
of biomolecular dynamics.  Before you start you should collect two .pdb files
you want to compare and move them into the DROIDS folder.  Naming convention
should be PDB_ID.pdb. Be sure to check that they are 'sensibly' homologous in
that they differ only in with regards to the effect you want to observe (i.e.
sequence difference, solvent or binding state).  Edit in Chimera if neccessary.
Remove crystalographic waters, mirrored structures or unusual ligands. Atypical
Amber preparations (e.g. beyond adding H and missing atoms using teleap) can be
done at the command line using Antechamber for further ligand library prep

NOTE: this program assumes .pdb files are completely ready for teLeap

Dependencies - perl, perl module (Descriptive), perl-tk, python, python-tk,
  python-gi, R-base, R-dev, R package(ggplot2), USCF Chimera 1.11, evince(pdf viewer)
  Amber16 (licensed from Univ of Ca; visit ambermd.org), Ambertools16
 (tested on Linux Mint 18.1 Cinnamon 64-bit Kernel 4.4.0-53-generic)

BabbittLab - Rochester Inst. Technol. Rochester NY

DROIDS 1.0               Copyright 2017 G.A. Babbitt.\n\n";


print "continue to GPL license? (y/n)\n";
my $go = <STDIN>;
chop($go);

if ($go eq "n") {exit;}

print " 

    DROIDS 1.0 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DROIDS 1.0 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DROIDS 1.0.  If not, see <http://www.gnu.org/licenses/>.

    Visit us on GitHub. \n\n";

print "continue to GUI ? (y/n)\n";
my $go = <STDIN>;
chop($go);

if ($go eq "n") {exit;}

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### Declare variables ####
my $fileIDq = '';
my $fileIDr = '';
my $forceID = '';
my $runsID = '';
my $implicit=0;
my $explicit=0;
my $solvType = '';
my $cutoffValueHeat=100;
my $cutoffValueEq=10;
my $cutoffValueProd=10;
my $cutoffValueSalt=0.0;
my $cutoffValueHeatFS=0;
my $cutoffValueEqFS=0;
my $cutoffValueProdFS=0;

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("AMBER MD control settings"); # Titles the main window

my $MDheatScale = $mw->Scale(-label=>"Length of MD heating run (ps) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>1000,
			-variable=>\$cutoffValueHeat,
			-tickinterval=>200,
			-resolution=>10,
			-length=>205
			);

my $MDeqScale = $mw->Scale(-label=>"Length of MD equilibration run (ns) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>100,
			-variable=>\$cutoffValueEq,
			-tickinterval=>20,
			-resolution=>1,
			-length=>205
			);

my $MDprodScale = $mw->Scale(-label=>"Length of each MD sample run (ps) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>1000,
			-variable=>\$cutoffValueProd,
			-tickinterval=>200,
			-resolution=>10,
			-length=>205
			);

my $MDsaltScale = $mw->Scale(-label=>"extra salt conc (M) (implicit only)  :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>0.6,
			-variable=>\$cutoffValueSalt,
			-tickinterval=>0.2,
			-resolution=>0.05,
			-length=>205
			);

# Solvation Frame
my $solnFrame = $mw->Frame(	-label => "METHOD OF SOLVATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $implicitCheck = $solnFrame->Radiobutton( -text => "implicit - Generalized Born",
						-value=>"im",
						-variable=>\$solvType
						);
	my $explicitCheck = $solnFrame->Radiobutton( -text => "explicit - Particle Mesh Ewald",
						-value=>"ex",
						-variable=>\$solvType
						);


# PDB ID Frame				
my $pdbFrame = $mw->Frame();
	my $QfileFrame = $pdbFrame->Frame();
		my $QfileLabel = $QfileFrame->Label(-text=>"pdb ID query (e.g. 1ubq) : ");
		my $QfileEntry = $QfileFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$fileIDq
					);
	my $RfileFrame = $pdbFrame->Frame();
		my $RfileLabel = $RfileFrame->Label(-text=>"pdb ID reference (e.g. 2ubq) : ");
		my $RfileEntry = $RfileFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$fileIDr
					);
		
	my $forceFrame = $pdbFrame->Frame();
		my $forceLabel = $forceFrame->Label(-text=>"Force Field (e.g. protein.ff14SB): ");
		my $forceEntry = $forceFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$forceID
					);
	my $runsFrame = $pdbFrame->Frame();
		my $runsLabel = $runsFrame->Label(-text=>"number of repeated MD sample runs: ");
		my $runsEntry = $runsFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$runsID
					);
		
# Buttons
my $controlButton = $mw -> Button(-text => "make MD control files (.ctl)", 
				-command => \&control
				); # Creates a ctl file button

my $launchButton = $mw -> Button(-text => "launch MD run (pmemd.cuda)", 
				-command => \&launch,
				-background => 'gray45',
				-foreground => 'white'
				); # Creates a launch button

my $killButton = $mw -> Button(-text => "kill MD run (pmemd.cuda) on GPU", 
				-command => \&kill
				); # Creates a kill button

my $survButton = $mw -> Button(-text => "open GPU job survellience", 
				-command => \&surv
				); # Creates a surv button

my $skipButton = $mw -> Button(-text => "go directly to CPPTRAJ (MD done previously)", 
				-command => \&skip
				); # Creates a skip button
my $teLeapButton = $mw -> Button(-text => "generate topology and coordinate files (teLeap)", 
				-command => \&teLeap
				); # Creates a teLeap button
my $alignButton = $mw -> Button(-text => "create sequence and structural alignment (UCSF Chimera)", 
				-command => \&align
				); # Creates a align button


#### Organize GUI Layout ####

$launchButton->pack(-side=>"right",
			-anchor=>"s"
			);
$teLeapButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$alignButton->pack(-side=>"bottom",
			-anchor=>"s"
    		);
$skipButton->pack(-side=>"bottom",
			-anchor=>"s"
    		);
$killButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$survButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$controlButton->pack(-side=>"bottom",
			-anchor=>"s"
			);

$QfileLabel->pack(-side=>"left");
$QfileEntry->pack(-side=>"left");
$RfileLabel->pack(-side=>"left");
$RfileEntry->pack(-side=>"left");
$forceLabel->pack(-side=>"left");
$forceEntry->pack(-side=>"left");
$runsLabel->pack(-side=>"left");
$runsEntry->pack(-side=>"left");

$forceFrame->pack(-side=>"top",
		-anchor=>"e");
$QfileFrame->pack(-side=>"top",
		-anchor=>"e");
$RfileFrame->pack(-side=>"top",
		-anchor=>"e");
$runsFrame->pack(-side=>"top",
		-anchor=>"e");
$pdbFrame->pack(-side=>"top",
		-anchor=>"n");

$implicitCheck->pack();
$explicitCheck->pack();
$solnFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$MDheatScale->pack(-side=>"top");
$MDeqScale->pack(-side=>"top");
$MDprodScale->pack(-side=>"top");
$MDsaltScale->pack(-side=>"top");

MainLoop; # Allows Window to Pop Up


########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################


sub control { # Write a control file and then call appropriate scripts that reference control file
	if ($solvType eq "im") {$repStr = "implicit";}
	if ($solvType eq "ex") {$repStr = "explicit";}
	
	# convert all times to femtosec
	$cutoffValueHeatFS = $cutoffValueHeat*1000;
	$cutoffValueEqFS = $cutoffValueEq*1000000;
	$cutoffValueProdFS = $cutoffValueProd*1000;

### make qury protein control file ###	
open(my $ctlFile1, '>', "MDq.ctl") or die "Could not open output file";
print $ctlFile1 
"PDB_ID\t$fileIDq\t# Protein Data Bank ID for MD run
Force_Field\t$forceID\t# AMBER force field to use in MD runs
Number_Runs\t$runsID\t# number of repeated samples of MD runs
Heating_Time\t$cutoffValueHeatFS\t# length of heating run (fs)
Equilibration_Time\t$cutoffValueEqFS\t# length of equilibration run (fs)
Production_Time\t$cutoffValueProdFS\t# length of production run (fs)
Solvation_Method\t$repStr\t# method of solvation (implicit or explicit)
Salt_Conc\t$cutoffValueSalt\t# salt concentration (implicit only, PME=O)";
close $ctlFile1;
### make qury protein control file ###	
open(my $ctlFile2, '>', "MDr.ctl") or die "Could not open output file";
print $ctlFile2 
"PDB_ID\t$fileIDr\t# Protein Data Bank ID for MD run
Force_Field\t$forceID\t# AMBER force field to use in MD runs
Number_Runs\t$runsID\t# number of repeated samples of MD runs
Heating_Time\t$cutoffValueHeatFS\t# length of heating run (fs)
Equilibration_Time\t$cutoffValueEqFS\t# length of equilibration run (fs)
Production_Time\t$cutoffValueProdFS\t# length of production run (fs)
Solvation_Method\t$repStr\t# method of solvation (implicit or explicit)
Salt_Conc\t$cutoffValueSalt\t# salt concentration (implicit only, PME=O)";
close $ctlFile2;

print "control files are made (see MDq.ctl and MDr.ctl)\n";

}


#####################################################################################################

sub teLeap { # create topology and coordinate files 
system "perl teLeap_proteinQuery.pl\n";
system "perl teLeap_proteinReference.pl\n";
}


######################################################################################################

sub launch { # launch MD run
system "perl MD_proteinQuery.pl\n";
system "perl MD_proteinReference.pl\n";
}

######################################################################################################

sub kill { # kill MD run
system "pkill pmemd\n";	
}

######################################################################################################

sub surv {
	### open job survalience terminals ######
system "x-terminal-emulator -e top\n";
system "x-terminal-emulator -e nvidia-smi -l 20\n";
}

######################################################################################################

sub skip{ # skip MD if .nc files are already generated
print "  skipping MD simulation and going straight to vector trajectory analysis\n\n";	
sleep (1);
system "perl GUI2_DROIDS.pl";
}

######################################################################################################

sub align{

print "STEP 1 - Here you will need to run MatchMaker in UCSF Chimera\n\n";
print "STEP 2 - Then run Match-Align in UCSF Chimera\n\n";
print "            if satisfied with alignment, save as a clustal file with ref PDB ID\n";
print "            in title (e.g. 1ubq_align.aln)\n\n";

print "continue? (y/n)\n";
my $go = <STDIN>;
chop($go);
if ($go eq "n") {exit;}
sleep(1);
print "            opening USCF Chimera and loading PDB ref structure\n\n";
print "            CREATE YOUR STRUCTURAL/SEQUENCE ALIGNMENT (.aln) NOW \n\n";
system("$chimera_path"."chimera $fileIDr.pdb $fileIDq.pdb\n");

sleep(0.5);
print "\n\n alignment procedure is complete\n";
sleep(0.5);
	
}


