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

    DROIDS 2.0 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DROIDS 2.0 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DROIDS 2.0.  If not, see <http://www.gnu.org/licenses/>.

    Visit us on GitHub. \n\n";

print "continue to GUI ? (y/n)\n";
my $go = <STDIN>;
chop($go);

if ($go eq "n") {exit;}

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- analytical engine and visual toolbox for functional evolutionary
  comparison of molecular dynamic simulation \n\n";

#### open PATHS specification ####

system "perl PATHS.pl\n";

#### Declare variables ####
my $testType = '';
my $gpuType = '';
my $mdType = '';

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("DROIDS 2.0 - free software for comparative protein dynamics"); # Titles the main window
$mw->setPalette("gray");

# Analysis Pipeline Frame
my $pipeFrame = $mw->Frame(	-label => "CHOOSE MODE OF MOLECULAR DYNAMIC COMPARISON",
				-relief => "groove",
				-borderwidth => 2
				);
	my $spRadio = $pipeFrame->Radiobutton( -text => "analyze impact of mutation(s) on protein or protein complexes (requires 2 PDB ID's with and w/o mutation)",
						-foreground => 'navy',
                        -value=>"sp",
						-variable=>\$testType
						);
    my $dsRadio = $pipeFrame->Radiobutton( -text => "diagnose impact of mutation in an unknown variant   (requires 4 PDB ID's = disease/unknown variant + human/ortholog pair)",
						-foreground => 'navy',
                        -value=>"ds",
						-variable=>\$testType
						);
    my $mpRadio = $pipeFrame->Radiobutton( -text => "analyze the rank of individual impacts of mutations on a protein  (requires 2 PDB ID's with and w/o mutation)",
						-foreground => 'navy',
                        -value=>"mp",
						-variable=>\$testType
						);
	my $dpRadio = $pipeFrame->Radiobutton( -text => "analyze impact of mutation(s) on DNA-protein interaction   (e.g. requires PDB ID's with and w/o cis or trans regulatory mutation)",
						-foreground => 'navy',
                        -value=>"dp",
						-variable=>\$testType
						);
    my $lpRadio = $pipeFrame->Radiobutton( -text => "analyze impact of a ligand interaction with a protein   (e.g. requires PDB ID's with and w/o drug, toxin, or activator)",
						-foreground => 'navy',
                        -value=>"lp",
						-variable=>\$testType
						);
	my $epRadio = $pipeFrame->Radiobutton( -text => "analyze impact of an epigenetic modification on a protein   (e.g requires PDB ID's with and w/o S-S bond or phosphorylation)",
						-foreground => 'navy',
                        -value=>"ep",
						-variable=>\$testType
						);
# System Type Frame
my $gpuFrame = $mw->Frame(	-label => "HOW MANY EFFECTIVE GPU's IN SYSTEM?",
				-relief => "groove",
				-borderwidth => 2
				);
	my $gpu1Radio = $gpuFrame->Radiobutton( -text => "single GPU system - run MD sequentially (single or multiple GPU with SLI)",
						-foreground => 'darkred',
                        -value=>"gpu1",
						-variable=>\$gpuType
						);
    my $gpu2Radio = $gpuFrame->Radiobutton( -text => "double GPU system - run query/ref MD simultaneously (i.e. no SLI)",
						-foreground => 'darkred',
                        -value=>"gpu2",
						-variable=>\$gpuType
						);

# MD Type Frame
my $mdFrame = $mw->Frame(	-label => "PRE-INSTALLED SOFTWARE FOR MD SIMULATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $md1Radio = $mdFrame->Radiobutton( -text => "AMBER 16 - licensed by Univ of CA",
						-foreground => 'darkgreen',
                        -value=>"amber16",
						-variable=>\$mdType
						);
    my $md2Radio = $mdFrame->Radiobutton( -text => "AMBER 18 - licensed by Univ of CA",
						-foreground => 'darkgreen',
                        -value=>"amber18",
						-variable=>\$mdType
						);
    my $md3Radio = $mdFrame->Radiobutton( -text => "OpenMM - open source MD sim library",
						-foreground => 'darkgreen',
                        -value=>"openMM",
						-variable=>\$mdType
						);

		
		
# Buttons

my $pipeButton = $mw -> Button(-text => "run DROIDS", 
				-command => \&go,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a go button
my $exitButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a go button



#### Organize GUI Layout ####
$pipeFrame->pack(-side=>"top",
		-anchor=>"s");

$spRadio->pack(-anchor=>"w");
$dsRadio->pack(-anchor=>"w");
$mpRadio->pack(-anchor=>"w");
$dpRadio->pack(-anchor=>"w");
$lpRadio->pack(-anchor=>"w");
$epRadio->pack(-anchor=>"w");

$gpuFrame->pack(-side=>"top",
		-anchor=>"s");
$gpu1Radio->pack(-anchor=>"w");
$gpu2Radio->pack(-anchor=>"w");

$mdFrame->pack(-side=>"top",
		-anchor=>"s");
$md1Radio->pack(-anchor=>"w");
$md2Radio->pack(-anchor=>"w");
$md3Radio->pack(-anchor=>"w");

$pipeButton->pack(-side=>"top",
			-anchor=>"s");
$exitButton->pack(-side=>"top",
			-anchor=>"s");

MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

sub stop {exit;}

########################################################################################
sub go {
# make control file for DROIDS	
print("launching DROIDS 2.0...\n");
  if ($mdType eq "amber16"){  
    if ($testType eq "sp" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSsp.pl\n";}  
	elsif ($testType eq "sp" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSsp_dualGPU.pl\n";}   
	elsif ($testType eq "ds" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ds" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "dp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "lp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "lp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSep.pl\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSep_dualGPU.pl\n";}
    else {print " PLEASE SELECT OPTIONS\n"}
  }
  
  if ($mdType eq "amber18"){  
    if ($testType eq "sp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}  
	elsif ($testType eq "sp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}   
	elsif ($testType eq "ds" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ds" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "dp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "lp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "lp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    else {print " PLEASE SELECT OPTIONS\n"}
  }
  
  if ($mdType eq "openMM"){  
    if ($testType eq "sp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}  
	elsif ($testType eq "sp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}   
	elsif ($testType eq "ds" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ds" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "mp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "dp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "lp" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "lp" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    elsif ($testType eq "ep" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
    else {print " PLEASE SELECT OPTIONS\n"}
  }
  
}

########################################################################
