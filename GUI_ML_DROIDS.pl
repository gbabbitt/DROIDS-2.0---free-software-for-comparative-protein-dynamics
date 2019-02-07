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
my $queryID = '';
my $refID = '';
my $lengthID = '';
my $homology = '';
my $chainN = '';


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
      #if ($header eq "start"){$startN = $value;}
      if ($header eq "homology"){$homology = $value;}
      if ($header eq "num_chains"){$chainN = $value;}
}
close IN;


#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("DROIDS - Machine Learning Classification"); # Titles the main window
$mw->setPalette("gray");

# ML method Frame
my $methodFrame = $mw->Frame(	-label => "MACHINE LEARNING METHOD",
				-relief => "groove",
				-borderwidth => 2
				);
	my $svmRadio = $methodFrame->Radiobutton( -text => "support vector machine",
						-foreground => 'navy',
                              -value=>"svm",
						-variable=>\$method
						);
	my $nnetRadio = $methodFrame->Radiobutton( -text => "neural network",
						-foreground => 'navy',
                              -value=>"nnet",
						-variable=>\$method
						);
	my $rforestRadio = $methodFrame->Radiobutton( -text => "random forest",
						-foreground => 'navy',
                              -value=>"rforest",
						-variable=>\$method
						);
# ML option Frame
my $optionFrame = $mw->Frame(	-label => "ENSEMBLE OPTIONS",
				-relief => "groove",
				-borderwidth => 2
				);
	my $noneRadio = $optionFrame->Radiobutton( -text => "none - single learner",
						-foreground => 'darkgreen',
                              -value=>"none",
						-variable=>\$option
						);
	my $bagRadio = $optionFrame->Radiobutton( -text => "bootstrap aggregation",
						-foreground => 'darkgreen',
                              -value=>"bag",
						-variable=>\$option
						);
	my $boostRadio = $optionFrame->Radiobutton( -text => "boosting",
						-foreground => 'darkgreen',
                              -value=>"boost",
						-variable=>\$option
						);
# movie viewing option Frame
my $viewFrame = $mw->Frame(	-label => "MOVIE VIEWING OPTIONS",
				-relief => "groove",
				-borderwidth => 2
				);
	my $monoRadio = $viewFrame->Radiobutton( -text => "monoscopic - normal 2D",
						-foreground => 'darkred',
                              -value=>"mono",
						-variable=>\$view
						);
	my $stereoRadio = $viewFrame->Radiobutton( -text => "stereoscopic - need 3D glasses",
						-foreground => 'darkred',
                              -value=>"stereo",
						-variable=>\$view
						);
	
# Buttons

my $stopButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a exit button
my $ctlButton = $mw -> Button(-text => "data preprocessing and control file", 
				-command => \&ctl,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $trainButton = $mw -> Button(-text => "training learner(s) on previous MD", 
				-command => \&train,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $mdButton = $mw -> Button(-text => "run new MD simulation to classify", 
				-command => \&md,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $deployButton = $mw -> Button(-text => "deploy learner(s) on new MD run", 
				-command => \&deploy,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $statsButton = $mw -> Button(-text => "return classification statistics", 
				-command => \&stats,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $movieButton = $mw -> Button(-text => "render dynamic classification movies", 
				-command => \&movie,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $playButton = $mw -> Button(-text => "play movies", 
				-command => \&play,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button

#### Organize GUI Layout ####
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$methodFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$optionFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$ctlButton->pack(-side=>"top",
			-anchor=>"s"
			);
$trainButton->pack(-side=>"top",
			-anchor=>"s"
			);
$mdButton->pack(-side=>"top",
			-anchor=>"s"
			);
$deployButton->pack(-side=>"top",
			-anchor=>"s"
			);
$statsButton->pack(-side=>"top",
			-anchor=>"s"
			);
$movieButton->pack(-side=>"top",
			-anchor=>"s"
			);
$viewFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$playButton->pack(-side=>"top",
			-anchor=>"s"
			);
$svmRadio->pack();
$nnetRadio->pack();
$rforestRadio->pack();
$noneRadio->pack();
$bagRadio->pack();
$boostRadio->pack();
$monoRadio->pack();
$stereoRadio->pack();
MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################
sub stop {exit;}
sub ctl {
     
}
sub train {
     
}
sub stats {
     
}
sub md {
     
}
sub deploy {
     
}
sub movie {
     
}
sub play {
if($view eq "mono") {print "\n\nquery state of protein is on left, reference state is on right\n\n";} 
if($view eq "stereo") {print "\n\nquery state of protein is on top, reference state is on bottom\n\n";} 
}







