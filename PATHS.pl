#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### Declare variables ####
my $chimera_path = '';
my $amber_path= '';
my $teleap_path = '';
my $openmm_path = '';

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("PATHS TO SOFTWARE (can exit if this ctl file already exists)"); # Titles the main window
$mw->setPalette("gray");


# PATH Frame				
my $pathFrame = $mw->Frame();
	my $amberFrame = $pathFrame->Frame();
		my $amberLabel = $amberFrame->Label(-text=>"path to amber home folder(e.g. /home/greg/Desktop/amber16/) : ");
		my $amberEntry = $amberFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$amber_path
					);
	my $chimeraFrame = $pathFrame->Frame();
		my $chimeraLabel = $chimeraFrame->Label(-text=>"path to chimera executable (e.g. /opt/UCSF/Chimera64-1.11/bin/) : ");
		my $chimeraEntry = $chimeraFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$chimera_path
					);
	my $teleapFrame = $pathFrame->Frame();
		my $teleapLabel = $teleapFrame->Label(-text=>"path to force fields (e.g. /home/greg/Desktop/amber16/dat/leap/cmd/) : ");
		my $teleapEntry = $teleapFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$teleap_path
					);	
	my $openmmFrame = $pathFrame->Frame();
		my $openmmLabel = $openmmFrame->Label(-text=>"path to OpenMM (e.g. NA ) : ");
		my $openmmEntry = $openmmFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$openmm_path
					);	
		
# Buttons
my $controlButton = $mw -> Button(-text => "make PATHS control file (.ctl)", 
				-command => \&control,
				-background => 'gray45',
                -foreground => 'white'
				); # Creates a ctl file button

my $exitButton = $mw -> Button(-text => "exit after PATHS are specified", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a go button

#### Organize GUI Layout ####
$exitButton->pack(-side=>"bottom",
			-anchor=>"s");
$controlButton->pack(-side=>"bottom",
			-anchor=>"s");

$amberLabel->pack(-side=>"left");
$amberEntry->pack(-side=>"left");
$chimeraLabel->pack(-side=>"left");
$chimeraEntry->pack(-side=>"left");
$teleapLabel->pack(-side=>"left");
$teleapEntry->pack(-side=>"left");
$openmmLabel->pack(-side=>"left");
$openmmEntry->pack(-side=>"left");

$amberFrame->pack(-side=>"top",
		-anchor=>"e");
$chimeraFrame->pack(-side=>"top",
		-anchor=>"e");
$teleapFrame->pack(-side=>"top",
		-anchor=>"e");
$openmmFrame->pack(-side=>"top",
		-anchor=>"e");
$pathFrame->pack(-side=>"top",
		-anchor=>"n");

MainLoop; # Allows Window to Pop Up


########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################
sub stop {exit;}
########################################################################################
sub control { # Write a control file and then call appropriate scripts that reference control file

### make qury protein control file ###	
open(my $ctlFile, '>', "paths.ctl") or die "Could not open output file";
print $ctlFile 
"amber_path\t$amber_path\t# path to amber home folder
chimera_path\t$chimera_path\t# path to Chimera executable
teleap_path\t$teleap_path\t# path to teLeap force field folder
openmm_path\t$openmm_path\t# path to OpenMM executable";
close $ctlFile;

print "control file for PATHS is done (see paths.ctl)\n";

}


#################################################################################




