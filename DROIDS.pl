#!/usr/bin/perl
use Tk;
use Tk::Text;
use Tk::Font;
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
        
   my $ssRadio = $pipeFrame->Radiobutton( -text => "(1) analyze stability (self-similarity) of dynamics on a single protein or protein complex									            (requires 1 PDB ID)",
						-foreground => 'navy',
                        -value=>"ss",
						-variable=>\$testType
						);
	 my $sdmRadio = $pipeFrame->Radiobutton( -text => "(2) analyze impact of one or several mutations (site-directed mutagenesis) on a protein or complex			            (requires 1 PDB ID and list of (locations/types) of mutation)",
						-foreground => 'navy',
                        -value=>"sdm",
						-variable=>\$testType
						);
	 my $edRadio = $pipeFrame->Radiobutton( -text => "(3) analyze impact of functional / evolutionary divergence on a protein or protein complex	  			        (requires 2 PDB ID's representing an paralog or ortholog pair)",
						-foreground => 'navy',
                        -value=>"ed",
						-variable=>\$testType
						);
    my $dp1Radio = $pipeFrame->Radiobutton( -text => "(4) analyze impact of DNA-protein interaction upon binding    					   	                      (requires 2 PDB ID's for DNA-protein complex and protein-only)",
						-foreground => 'navy',
                        -value=>"dp1",
						-variable=>\$testType
						);
	 my $dp2Radio = $pipeFrame->Radiobutton( -text => "(5) analyze impact of mutation(s) on DNA-protein interaction		 		                 (requires 2 PDB ID's as in #4 and (locations/types) of cis or trans regulatory mutation)",
						-foreground => 'navy',
                        -value=>"dp2",
						-variable=>\$testType
						);
	 my $dp3Radio = $pipeFrame->Radiobutton( -text => "(6) analyze comparison of two DNA-protein interactions 				 			    	             (requires 2 PDB ID representing binding protein homologs)",
						-foreground => 'navy',
                        -value=>"dp3",
						-variable=>\$testType
						);
    my $lp1Radio = $pipeFrame->Radiobutton( -text => "(7)  analyze impact of a protein-ligand interaction upon binding drug, toxin, or activator  				(requires 3 PDB ID's = protein-ligand complex, protein-only, and ligand only)",
						-foreground => 'navy',
                        -value=>"lp1",
						-variable=>\$testType
						);
	 my $lp2Radio = $pipeFrame->Radiobutton( -text => "(8)  analyze impact of mutation(s) on protein-ligand interaction		 			   	            (requires 3 PDB ID's as in #7 and list of (locations/types) of mutation)",
						-foreground => 'navy',
                        -value=>"lp2",
						-variable=>\$testType
						);
	 
# System Type Frame
my $gpuFrame = $mw->Frame(	-label => "HOW MANY EFFECTIVE GPU's IN SYSTEM?",
				-relief => "groove",
				-borderwidth => 2
				);
	my $gpu1Radio = $gpuFrame->Radiobutton( -text => "single GPU system - run MD sequentially (single GPU)                           NOTE: can only run single chain protein",
						-foreground => 'darkred',
                        -value=>"gpu1",
						-variable=>\$gpuType
						);
        my $gpu2Radio = $gpuFrame->Radiobutton( -text => "double GPU system - run query/ref MD simultaneously (i.e. no SLI)    NOTE: can also analyze multi-chain proteins ",
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

$ssRadio->pack(-anchor=>"w");
$sdmRadio->pack(-anchor=>"w");
$edRadio->pack(-anchor=>"w");

$dp1Radio->pack(-anchor=>"w");
$dp2Radio->pack(-anchor=>"w");
$dp3Radio->pack(-anchor=>"w");
$lp1Radio->pack(-anchor=>"w");
$lp2Radio->pack(-anchor=>"w");


$gpuFrame->pack(-side=>"top",
		-anchor=>"s");
$gpu1Radio->pack(-anchor=>"w");
$gpu2Radio->pack(-anchor=>"w");

$mdFrame->pack(-side=>"top",
		-anchor=>"s");
$md1Radio->pack(-anchor=>"w");
$md2Radio->pack(-anchor=>"w");
#$md3Radio->pack(-anchor=>"w");

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
     if ($testType eq "ss" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSss.pl\n";}  
	elsif ($testType eq "ss" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSss_dualGPU.pl\n";}
   elsif ($testType eq "sdm" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSsdm.pl\n";}  
	elsif ($testType eq "sdm" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSsdm_dualGPU.pl\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSed.pl\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSed_dualGPU.pl\n";}
	elsif ($testType eq "dp1" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSdp1.pl\n";}
        elsif ($testType eq "dp1" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSdp1_dualGPU.pl\n";}
	elsif ($testType eq "dp2" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSdp2.pl\n";}
        elsif ($testType eq "dp2" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSdp2_dualGPU.pl\n";}
	elsif ($testType eq "dp3" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSdp3.pl\n";}
        elsif ($testType eq "dp3" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSdp3_dualGPU.pl\n";}
	elsif ($testType eq "lp1" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSlp1.pl\n";}
        elsif ($testType eq "lp1" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSlp1_dualGPU.pl\n";}
	elsif ($testType eq "lp2" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSlp2.pl\n";}
        elsif ($testType eq "lp2" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSlp2_dualGPU.pl\n";}
        
        else {print " PLEASE SELECT OPTIONS\n"}
        }
  
  if ($mdType eq "amber18"){  
     if ($testType eq "ss" && $gpuType eq "gpu1") {print " AMBER 18 will be supported soon\n";}  
	elsif ($testType eq "ss" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "sdm" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}  
	elsif ($testType eq "sdm" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu2") {print"  AMBER 18 will be supported soon\n";}
		elsif ($testType eq "dp1" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
        elsif ($testType eq "dp1" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "dp2" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
        elsif ($testType eq "dp2" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "dp3" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
        elsif ($testType eq "dp3" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
		  elsif ($testType eq "lp1" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
        elsif ($testType eq "lp1" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
	elsif ($testType eq "lp2" && $gpuType eq "gpu1") {print "  AMBER 18 will be supported soon\n";}
        elsif ($testType eq "lp2" && $gpuType eq "gpu2") {print "  AMBER 18 will be supported soon\n";}
        
        else {print " PLEASE SELECT OPTIONS\n"}
        }
  
  if ($mdType eq "openMM"){  
     if ($testType eq "ss" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}  
	elsif ($testType eq "ss" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "sdm" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}  
	elsif ($testType eq "sdm" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "ed" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "ds1" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "ds1" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp1" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "dp1" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp2" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "dp2" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "dp3" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "dp3" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
		  elsif ($testType eq "lp1" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "lp1" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
	elsif ($testType eq "lp2" && $gpuType eq "gpu1") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        elsif ($testType eq "lp2" && $gpuType eq "gpu2") {print " THIS OPTION IS NOT YET AVAILABLE\n";}
        
        else {print " PLEASE SELECT OPTIONS\n"}
        }
  
}

########################################################################
