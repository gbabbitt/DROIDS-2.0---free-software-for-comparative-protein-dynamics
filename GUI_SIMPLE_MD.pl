#!/usr/bin/perl
use Tk;
use strict;
use warnings;
use feature ":5.10";

#### Introductory Message #############

print "Welcome to simple MD script of DROIDS - this is intended for control
of GPU accelerated AMBER16 when doing MD on a single protein.  Use GUI_DROIDS.pl
if you intend to run a pair of simulations and compare them.  

Dependency - Perl, Perl/TK, Python, R, USCF Chimera, Amber16, Ambertools16
                          (tested on Linux Mint)

BabbittLab - Rochester Inst. Technol. \n\n";

print "continue? (y/n)\n";
my $go = <STDIN>;
chop($go);

if ($go eq "n") {exit;}


#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### Declare variables ####
my @rep;
my $fileID = '';
my $forceID = '';
my $runsID = '';
my $implicit=0;
my $explicit=0;
my $cutoffValueHeat=100;
my $cutoffValueEq=10;
my $cutoffValueProd=10;
my $cutoffValueHeatFS=0;
my $cutoffValueEqFS=0;
my $cutoffValueProdFS=0;
my $molType="dna";

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("AMBER MD control settings"); # Titles the main window

my $MDheatScale = $mw->Scale(-label=>"Length of MD heating (ps) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>1000,
			-variable=>\$cutoffValueHeat,
			-tickinterval=>200,
			-resolution=>10,
			-length=>200
			);

my $MDeqScale = $mw->Scale(-label=>"Length of MD equilibration (ns) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>100,
			-variable=>\$cutoffValueEq,
			-tickinterval=>20,
			-resolution=>1,
			-length=>200
			);

my $MDprodScale = $mw->Scale(-label=>"Length of each MD production (ns) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>100,
			-variable=>\$cutoffValueProd,
			-tickinterval=>20,
			-resolution=>1,
			-length=>200
			);

# Solvation Frame
my $solnFrame = $mw->Frame(	-label => "Method of Solvation",
				-relief => "groove",
				-borderwidth => 2
				);
	my $implicitCheck = $solnFrame->Checkbutton( -text => "implicit",
						-variable=>\$implicit
						);
	my $explicitCheck = $solnFrame->Checkbutton( -text => "explicit",
						-variable=>\$explicit
						);


# PDB ID Frame				
my $pdbFrame = $mw->Frame();
	my $fileFrame = $pdbFrame->Frame();
		my $fileLabel = $fileFrame->Label(-text=>"pdb ID (e.g. 5ire) : ");
		my $fileEntry = $fileFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$fileID
					);
	my $forceFrame = $pdbFrame->Frame();
		my $forceLabel = $forceFrame->Label(-text=>"ForceField (e.g. protein.ff14SB): ");
		my $forceEntry = $forceFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$forceID
					);
	my $runsFrame = $pdbFrame->Frame();
		my $runsLabel = $runsFrame->Label(-text=>"number of repeated MD runs: ");
		my $runsEntry = $runsFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$runsID
					);
		
# Buttons
my $controlButton = $mw -> Button(-text => "make MD.ctl file", 
				-command => \&control
				); # Creates a ctl file button

my $launchButton = $mw -> Button(-text => "launch MD run (pmemd.cuda)", 
				-command => \&launch
				); # Creates a launch button

my $killButton = $mw -> Button(-text => "kill MD run (pmemd.cuda) on GPU", 
				-command => \&kill
				); # Creates a kill button

my $survButton = $mw -> Button(-text => "open GPU job survellience", 
				-command => \&surv
				); # Creates a surv button


#### Organize GUI Layout ####

$launchButton->pack(-side=>"right",
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

$fileLabel->pack(-side=>"left");
$fileEntry->pack(-side=>"left");
$forceLabel->pack(-side=>"left");
$forceEntry->pack(-side=>"left");
$runsLabel->pack(-side=>"left");
$runsEntry->pack(-side=>"left");

$forceFrame->pack(-side=>"top",
		-anchor=>"e");
$fileFrame->pack(-side=>"top",
		-anchor=>"e");
$runsFrame->pack(-side=>"top",
		-anchor=>"e");
$pdbFrame->pack(-side=>"top",
		-anchor=>"n");


$MDheatScale->pack(-side=>"top");
$MDeqScale->pack(-side=>"top");
$MDprodScale->pack(-side=>"top");

$implicitCheck->pack();
$explicitCheck->pack();
$solnFrame->pack(-side=>"top",
		-anchor=>"n"
		);


MainLoop; # Allows Window to Pop Up


sub control { # Write a control file and then call appropriate scripts that reference control file
	if ($implicit == 1) {
		push @rep, "implicit";
	}
	if ($explicit == 1) {
		push @rep, "explicit";
	}
	my $repStr = join(",",@rep);
	
	# convert all times to femtosec
	$cutoffValueHeatFS = $cutoffValueHeat*1000;
	$cutoffValueEqFS = $cutoffValueEq*1000000;
	$cutoffValueProdFS = $cutoffValueProd*1000000;
	
	open(my $ctlFile, '>', "MD.ctl") or die "Could not open output file";

	print $ctlFile 
"PDB_ID\t$fileID\t# Protein Data Bank ID for MD run
Force_Field\t$forceID\t# AMBER force field to use in MD runs
Number_Runs\t$runsID\t# number of repeated samples of MD runs
Heating_Time\t$cutoffValueHeatFS\t# length of heating run (fs)
Equilibration_Time\t$cutoffValueEqFS\t# length of equilibration run (fs)
Production_Time\t$cutoffValueProdFS\t# length of production run (fs)
Solvation_Method\t$repStr\t# method of solvation (implicit or explicit)";

print "control file is made (see MD.ctl)\n";


}

sub launch { # launch MD run
system "perl MD_protein.pl\n";	
}

sub kill { # kill MD run
system "pkill pmemd\n";	

}

sub surv {
	### open job survalience terminals ######
system "x-terminal-emulator -e top\n";
system "x-terminal-emulator -e nvidia-smi -l 20\n";
}