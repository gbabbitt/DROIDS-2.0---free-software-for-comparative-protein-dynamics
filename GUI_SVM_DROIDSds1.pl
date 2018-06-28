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
#open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
#my @IN = <IN>;
#for (my $i = 0; $i < scalar @IN; $i++){
#	 my $INrow = $IN[$i];
#	 my @INrow = split (/\s+/, $INrow);
#	 my $header = @INrow[0];
#	 my $value = @INrow[1];
#	 if ($header eq "query"){$queryID = $value;}
#      if ($header eq "reference"){$refID = $value;}
#      if ($header eq "length"){$lengthID = $value;}
#      if ($header eq "start"){$startN = $value;}
#      if ($header eq "homology"){$homology = $value;}
#}
#close IN;


#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("SVM parameters"); # Titles the main window
$mw->setPalette("gray");

my $degreeScale = $mw->Scale(-label=>"degree (for polynomial) :",
			-orient=>'h',
			-digit=>3,
			-from=>1,
			-to=>5,
			-variable=>\$degree,
			-tickinterval=>1,
			-resolution=>1,
			-length=>200,
			);

my $gammaScale = $mw->Scale(-label=>"gamma (for radial function) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>12,
			-variable=>\$gamma,
			-tickinterval=>3,
			-resolution=>0.1,
			-length=>200,
			);


# SVM parameter Frame
my $svmFrame = $mw->Frame(	-label => "CHOOSE SVM KERNEL METHOD",
				-relief => "groove",
				-borderwidth => 2
				);
	my $linRadio = $svmFrame->Radiobutton( -text => "linear",
						-value=>"linear",
						-variable=>\$kernel
						);
	my $polyRadio = $svmFrame->Radiobutton( -text => "polynomial",
						-value=>"polynomial",
						-variable=>\$kernel
						);
	my $rdfRadio = $svmFrame->Radiobutton( -text => "radial basis function",
						-value=>"radial",
						-variable=>\$kernel
						);
	
# Buttons

my $ctlButton = $mw -> Button(-text => "write SVM control file (DROIDSsvm.ctl)", 
				-command => \&ctl,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $stopButton = $mw -> Button(-text => "exit back to stats analysis", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button

#### Organize GUI Layout ####


$svmFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$degreeScale->pack(-side=>"top");
$gammaScale->pack(-side=>"top");
$ctlButton->pack(-side=>"top",
			-anchor=>"s"
			);
$linRadio->pack();
$polyRadio->pack();
$rdfRadio->pack();
$stopButton->pack(-side=>"top",
			-anchor=>"s"
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

	$testStr = "flux"; $testStrLong = "fluctuation";  # file and folder labels
	
open(CTL, '>', "DROIDSsvm.ctl") or die "Could not open output file";
print CTL "kernel\t"."$kernel\t # kernel method \n";
print CTL "degree\t"."$degree\t # degree polynomial\n";
print CTL "gamma\t"."$gamma\t # gamma for rdf\n";
close CTL;
print("CTL file made\n");
}

#############################################################################################



