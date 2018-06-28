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
      if ($header eq "homology"){$homology = $value;}
}
close IN;


#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("DROIDS - Statistical Analysis"); # Titles the main window
$mw->setPalette("gray");

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
# local mutation defn Frame				
my $mutFrame = $mw->Frame();
	my $MfileFrame = $mutFrame->Frame();
		my $MsizeLabel = $mutFrame->Label(-text=>"define local mutation effect size for rank analysis\n (e.g. 5 = analyze dFLUX within 5 AA's of subs site): ");
		my $MsizeEntry = $mutFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$mut_targetsize
					);
		
# Buttons

my $statsButton = $mw -> Button(-text => "make statistical comparisons and plots in R for human/ortholog pair", 
				-command => \&stats,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a stats test button
my $stats2Button = $mw -> Button(-text => "make statistical comparisons and plots in R for human/disease variant pair", 
				-command => \&stats2,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a stats test button
my $stats3Button = $mw -> Button(-text => "make statistical comparisons and plots in R for human/novel variant pair", 
				-command => \&stats3,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a stats test button
my $nextButton = $mw -> Button(-text => "goto DROIDS image and movie rendering", 
				-command => \&next,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a next GUI button
my $stopButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a exit button
my $ctlButton = $mw -> Button(-text => "append control file (DROIDS.ctl)", 
				-command => \&ctl,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $rankButton = $mw -> Button(-text => "compare local impacts of disease vs orthologous mutations", 
				-command => \&rank,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $novelmdButton = $mw -> Button(-text => "run MD to classify impact of a set of novel mutations", 
				-command => \&novelMD,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $novelButton = $mw -> Button(-text => "run SVM to classify impact of the set of novel mutations", 
				-command => \&novelSVM,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button


#### Organize GUI Layout ####
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$nextButton->pack(-side=>"bottom",
			-anchor=>"s"
			);


$cutoffScale->pack(-side=>"top");
$mtcFrame->pack(-side=>"top",
		-anchor=>"s"
		);

$ctlButton->pack(-side=>"top",
			-anchor=>"s"
			);
$statsButton->pack(-side=>"top",
			-anchor=>"s"
			);
$stats2Button->pack(-side=>"top",
			-anchor=>"s"
			);
$mutFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$MsizeLabel->pack(-side=>"left");
$MsizeEntry->pack(-side=>"left");
$rankButton->pack(-side=>"top",
			-anchor=>"s"
			);
$novelmdButton->pack(-side=>"top",
			-anchor=>"s"
			);
$stats3Button->pack(-side=>"top",
			-anchor=>"s"
			);
$novelButton->pack(-side=>"top",
			-anchor=>"s"
			);
$noneRadio->pack();
$bonfRadio->pack();
$fdrRadio->pack();


MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################
sub stop {exit;}
########################################################################################
sub next {system "perl GUI_IMAGE_DROIDSds1.pl";}
########################################################################################
sub ctl {
# update control file for DROIDS	
print("updating ctl file...\n");

	$testStr = "flux"; $testStrLong = "fluctuation";  # file and folder labels
	
open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "test_type\t"."$testStr\t # test method\n";
close CTL;

open(CTL, '>>', "DROIDSmut.ctl") or die "Could not open output file";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "test_type\t"."$testStr\t # test method\n";
close CTL;

open(CTL, '>>', "DROIDSnovel.ctl") or die "Could not open output file";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "test_type\t"."$testStr\t # test method\n";
close CTL;

print("CTL file updated\n");


}

########################################################################

sub rank{
     
# reread control files
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
}
close IN;

##############
print "collecting local atom fluctuations near points of ortholog mutation\n";
sleep(1);
open(OUT1, ">"."mutate_flux_ortholog.txt") or die "could not open output file\n";
print OUT1 "mut_number\t"."subsLOCATION\t"."posAA\t"."refAA\t"."queryAA\t"."atomTYPE\t"."flux_ref\t"."flux_query\n";
open(MUT, "<"."mutate_list_ortholog.txt") or die "could not find list of mutations\n";
my @MUT = <MUT>;
for (my $m = 0; $m <= scalar @MUT; $m++){
  my $MUTrow = $MUT[$m];
  my @MUTrow = split(/\s+/, $MUTrow);
  my $subsTYPE = $MUTrow[0];
  if ($m => 1){$subsLOCATION = $MUTrow[1]}
  my $mut_number = $m;
  if ($subsLOCATION =~ m/\d/){print "mutation is near "."$subsLOCATION\n"};
  for (my $r = 0; $r <= $lengthID; $r++){
   
   $filenumber = $startN + $r;
   open(INFO, "<"."atomflux_ortholog/DROIDSfluctuation_$filenumber.txt") or next;
   my @INFO = <INFO>;
   for (my $rr = 0; $rr <= scalar @INFO; $rr++){
   my $INFOrow = $INFO[$rr];
	        my @INFOrow = split(/\s+/, $INFOrow); 
	        my $posAA = $INFOrow[1];
		   my $refAA = $INFOrow[2];
		   my $queryAA = $INFOrow[3];
             my $atomTYPE = $INFOrow[5];
             my $flux_ref = $INFOrow[6];
		   my $flux_query = $INFOrow[7];
     #print "$subsLOCATION\t"."$posAA\n";
     #if ($subLOCATION - $posAA <= $mut_targetsize){print "$mut_number\t"."$posAA\t"."$refAA\t"."$atomTYPE\t"."$flux_ref\t"."$$flux_query\n"}
     if ($subsLOCATION =~ m/\d/ && abs($subsLOCATION - $posAA) <=  $mut_targetsize){print OUT1 "$mut_number\t"."$subsLOCATION\t"."$posAA\t"."$refAA\t"."$queryAA\t"."$atomTYPE\t"."$flux_ref\t"."$flux_query\n"}
    }
    close INFO;   
  }
}
close MUT;
close OUT1;     
print "done collecting local atom fluctuations near points of ortholog mutation\n\n";
sleep(1);
print "analyzing local atom fluctuations near points of ortholog mutation\n";
sleep(1);
open(IN, "<"."DROIDS.ctl") or die "could not find CPPTRAJ input control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "query") { $queryID = $value;}
    if ($header eq "cutoff_value") { $level_sig = $value;}
    if ($header eq "test_type") { $seq_struct_flux = $value;}
close IN;
}

open(OUT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt")or die "could not open statistics.txt\n";
print OUT2 "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\n";

open(OUT3, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations_diseaseANDortholog.txt")or die "could not open statistics.txt\n";
print OUT3 "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\t"."class\n";

@fluxR = ();
@fluxQ = ();
open(MUT2, "<"."mutate_flux_ortholog.txt") or die "could not find list of mutations\n";
my @MUT2 = <MUT2>;
for (my $s = 0; $s <= scalar @MUT2; $s++){
  if ($s == 0){next;}
  my $MUTrow = $MUT2[$s];
  my $nextMUTrow = $MUT2[$s+1];
  my @MUTrow = split(/\s+/, $MUTrow);
  my @nextMUTrow = split(/\s+/, $nextMUTrow);
  my $mut_number = $MUTrow[0];
  my $next_number = $nextMUTrow[0];
  my $mut_location = $MUTrow[1];
  my $test_location = $MUTrow[2];
  my $ref_AA = $MUTrow[3];
  my $mut_AA = $MUTrow[4];
  my $mut_fluxR = $MUTrow[6];
  my $mut_fluxQ = $MUTrow[7];
  push (@fluxR, $mut_fluxR);
  push (@fluxQ, $mut_fluxQ);
  #print "$mut_number\t"."$next_number\n";
  if ($mut_location == $test_location){$label = "$ref_AA"."->"."$mut_AA"."$mut_location";}
  if ($mut_number != $next_number){
     print scalar @fluxR; print "\n";
     print scalar @fluxQ; print "\n";
     
     # calculate avg dFLUX
     $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
     $statSCORE->add_data (@fluxR);
	$flux_ref_avg = $statSCORE->mean();
     #$flux_ref_n = $statSCORE->count();
     #print "flux_ref_n\t"."$flux_ref_n\n";
	$statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
     $statSCORE->add_data (@fluxQ);
     $flux_query_avg = $statSCORE->mean();
     #$flux_query_n = $statSCORE->count();
     #print "flux_query_n\t"."$flux_query_n\n";
     $delta_flux = ($flux_ref_avg - $flux_query_avg);
     # calculate JS divergence
     open (TMP1, ">"."flux_mut_temp.txt") or die "could not create temp file\n";
     print TMP1 "flux_ref\t"."flux_query\n";
     for (my $t = 0; $t <= scalar @fluxQ; $t++){print TMP1 "$fluxR[$t]\t"; print TMP1 "$fluxQ[$t]\n";}
     close TMP1;
     open (TMP2, ">"."flux_mut_KL.txt") or die "could not create temp file\n";
     close TMP2;
     open (Rinput, "| R --vanilla")||die "could not start R command line\n";
     print Rinput "library('FNN')\n";
     print Rinput "data = read.table('flux_mut_temp.txt', header = TRUE)\n"; 
     $flux_ref = "data\$flux_ref"; # flux on reference residue
     $flux_query = "data\$flux_query"; # flux on query residue
     print Rinput "d1 = data.frame(fluxR=$flux_ref, fluxQ=$flux_query)\n";
     #print Rinput "print(d1)\n";
     print Rinput "myKL<-KL.dist($flux_ref, $flux_query, k=10)\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink('flux_mut_KL.txt')\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink()\n";
     # write to output file and quit R
     print Rinput "q()\n";# quit R 
     print Rinput "n\n";# save workspace image?
     close Rinput;
     open (TMP3, "<"."flux_mut_KL.txt") or die "could not create temp file\n";
     my @TMP3 = <TMP3>;
     for (my $tt = 0; $tt <= scalar @TMP3; $tt++){
     $TMP3row = $TMP3[$tt];
     @TMP3row = split (/\s+/, $TMP3row);
     $header = $TMP3row[0];
     $value = $TMP3row[1];
     #print "$header\t"."$value\n";
     if ($header eq "[1]"){$KL = $value;}
     }
     if ($delta_flux >= 0){$effect = "PROTEIN_STABILIZED";} # make KL value negative if dFLUX is negative
     if ($delta_flux < 0){$effect = "PROTEIN_DESTABILIZED";} # make KL value negative if dFLUX is negative
     print "my KL is "."$KL\n";
     close TMP3;
     print OUT2 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\n";
     print OUT3 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\t"."ortholog\n";
     @fluxR = ();
     @fluxQ = ();
     }
 
}
close MUT2;
close OUT2;
sleep(1);
print "\ndone analyzing local atom fluctuations near points of ortholog mutation\n\n";
sleep(1);
####################

# reset query to disease
open(IN, "<"."DROIDSmut.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      
}
close IN;


print "collecting local atom fluctuations near points of disease mutation\n";
sleep(1);
open(OUT1, ">"."mutate_flux_disease.txt") or die "could not open output file\n";
print OUT1 "mut_number\t"."subsLOCATION\t"."posAA\t"."refAA\t"."queryAA\t"."atomTYPE\t"."flux_ref\t"."flux_query\n";
open(MUT, "<"."mutate_list_disease.txt") or die "could not find list of mutations\n";
my @MUT = <MUT>;
for (my $m = 0; $m <= scalar @MUT; $m++){
  my $MUTrow = $MUT[$m];
  my @MUTrow = split(/\s+/, $MUTrow);
  my $subsTYPE = $MUTrow[0];
  if ($m => 1){$subsLOCATION = $MUTrow[1]}
  my $mut_number = $m;
  if ($subsLOCATION =~ m/\d/){print "mutation is near "."$subsLOCATION\n"};
  for (my $r = 0; $r <= $lengthID; $r++){
   
   $filenumber = $startN + $r;
   open(INFO, "<"."atomflux_disease/DROIDSfluctuation_$filenumber.txt") or next;
   my @INFO = <INFO>;
   for (my $rr = 0; $rr <= scalar @INFO; $rr++){
   my $INFOrow = $INFO[$rr];
	        my @INFOrow = split(/\s+/, $INFOrow); 
	        my $posAA = $INFOrow[1];
		   my $refAA = $INFOrow[2];
		   my $queryAA = $INFOrow[3];
             my $atomTYPE = $INFOrow[5];
             my $flux_ref = $INFOrow[6];
		   my $flux_query = $INFOrow[7];
     #print "$subsLOCATION\t"."$posAA\n";
     #if ($subLOCATION - $posAA <= $mut_targetsize){print "$mut_number\t"."$posAA\t"."$refAA\t"."$atomTYPE\t"."$flux_ref\t"."$$flux_query\n"}
     if ($subsLOCATION =~ m/\d/ && abs($subsLOCATION - $posAA) <=  $mut_targetsize){print OUT1 "$mut_number\t"."$subsLOCATION\t"."$posAA\t"."$refAA\t"."$queryAA\t"."$atomTYPE\t"."$flux_ref\t"."$flux_query\n"}
    }
    close INFO;   
  }
}
close MUT;
close OUT1;     
print "done collecting local atom fluctuations near points of disease mutation\n\n";
sleep(1);
print "analyzing local atom fluctuations near points of disease mutation\n";
sleep(1);
open(IN, "<"."DROIDS.ctl") or die "could not find CPPTRAJ input control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "query") { $queryID = $value;}
    if ($header eq "cutoff_value") { $level_sig = $value;}
    if ($header eq "test_type") { $seq_struct_flux = $value;}
close IN;
}

# reset query to disease
open(IN, "<"."DROIDSmut.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      
}
close IN;

open(OUT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt")or die "could not open statistics.txt\n";
print OUT2 "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\n";

@fluxR = ();
@fluxQ = ();
open(MUT2, "<"."mutate_flux_disease.txt") or die "could not find list of mutations\n";
my @MUT2 = <MUT2>;
for (my $s = 0; $s <= scalar @MUT2; $s++){
  if ($s == 0){next;}
  my $MUTrow = $MUT2[$s];
  my $nextMUTrow = $MUT2[$s+1];
  my @MUTrow = split(/\s+/, $MUTrow);
  my @nextMUTrow = split(/\s+/, $nextMUTrow);
  my $mut_number = $MUTrow[0];
  my $next_number = $nextMUTrow[0];
  my $mut_location = $MUTrow[1];
  my $test_location = $MUTrow[2];
  my $ref_AA = $MUTrow[3];
  my $mut_AA = $MUTrow[4];
  my $mut_fluxR = $MUTrow[6];
  my $mut_fluxQ = $MUTrow[7];
  push (@fluxR, $mut_fluxR);
  push (@fluxQ, $mut_fluxQ);
  #print "$mut_number\t"."$next_number\n";
  if ($mut_location == $test_location){$label = "$ref_AA"."->"."$mut_AA"."$mut_location";}
  if ($mut_number != $next_number){
     print scalar @fluxR; print "\n";
     print scalar @fluxQ; print "\n";
     
     # calculate avg dFLUX
     $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
     $statSCORE->add_data (@fluxR);
	$flux_ref_avg = $statSCORE->mean();
     #$flux_ref_n = $statSCORE->count();
     #print "flux_ref_n\t"."$flux_ref_n\n";
	$statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
     $statSCORE->add_data (@fluxQ);
     $flux_query_avg = $statSCORE->mean();
     #$flux_query_n = $statSCORE->count();
     #print "flux_query_n\t"."$flux_query_n\n";
     $delta_flux = ($flux_ref_avg - $flux_query_avg);
     # calculate JS divergence
     open (TMP1, ">"."flux_mut_temp.txt") or die "could not create temp file\n";
     print TMP1 "flux_ref\t"."flux_query\n";
     for (my $t = 0; $t <= scalar @fluxQ; $t++){print TMP1 "$fluxR[$t]\t"; print TMP1 "$fluxQ[$t]\n";}
     close TMP1;
     open (TMP2, ">"."flux_mut_KL.txt") or die "could not create temp file\n";
     close TMP2;
     open (Rinput, "| R --vanilla")||die "could not start R command line\n";
     print Rinput "library('FNN')\n";
     print Rinput "data = read.table('flux_mut_temp.txt', header = TRUE)\n"; 
     $flux_ref = "data\$flux_ref"; # flux on reference residue
     $flux_query = "data\$flux_query"; # flux on query residue
     print Rinput "d1 = data.frame(fluxR=$flux_ref, fluxQ=$flux_query)\n";
     #print Rinput "print(d1)\n";
     print Rinput "myKL<-KL.dist($flux_ref, $flux_query, k=10)\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink('flux_mut_KL.txt')\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink()\n";
     # write to output file and quit R
     print Rinput "q()\n";# quit R 
     print Rinput "n\n";# save workspace image?
     close Rinput;
     open (TMP3, "<"."flux_mut_KL.txt") or die "could not create temp file\n";
     my @TMP3 = <TMP3>;
     for (my $tt = 0; $tt <= scalar @TMP3; $tt++){
     $TMP3row = $TMP3[$tt];
     @TMP3row = split (/\s+/, $TMP3row);
     $header = $TMP3row[0];
     $value = $TMP3row[1];
     #print "$header\t"."$value\n";
     if ($header eq "[1]"){$KL = $value;}
     }
     if ($delta_flux >= 0){$effect = "PROTEIN_STABILIZED";} # make KL value negative if dFLUX is negative
     if ($delta_flux < 0){$effect = "PROTEIN_DESTABILIZED";} # make KL value negative if dFLUX is negative
     print "my KL is "."$KL\n";
     close TMP3;
     print OUT2 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\n";
     print OUT3 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\t"."disease\n";
     @fluxR = ();
     @fluxQ = ();
     }
 
}
close MUT2;
close OUT2;
close OUT3;
sleep(1);
print "\ndone analyzing local atom fluctuations near points of disease mutation\n\n";
sleep(1);

####################

# reset values from main control file
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      
}
close IN;

####################
print "\nplotting bar chart of local impacts near points of mutation\n\n";
sleep(1);
open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";

# read data into R
print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt', header = TRUE)\n"; 

$location1 = "data1\$location"; # position 	
$label1 = "data1\$label"; # AA label 
$effect1 = "data1\$effect"; # effect of mutation
$KLdivergence1 = "data1\$KLdivergence"; # impact of mutation
$queryAA1 = "data1\$startqueryAA"; # avg flux on query structure
print Rinput "data1\n";
print Rinput "data2 = read.table('./DROIDS_results_$refID"."mut"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt', header = TRUE)\n"; 
$location2 = "data2\$location"; # position 	
$label2 = "data2\$label"; # AA label 
$effect2 = "data2\$effect"; # effect of mutation
$KLdivergence2 = "data2\$KLdivergence"; # impact of mutation
$queryAA2 = "data2\$startqueryAA"; # avg flux on query structure
print Rinput "data2\n";
########## global KS test #### does disease mutation significantly differ from orthologs
print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
print Rinput "print(ks_test)\n";
#print to file
open(STAT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStest_diseaseVSortholog.txt")or die "could not open statistics.txt\n";
close STAT;
print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStest_diseaseVSortholog.txt')\n";
print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
print Rinput "print(ks_test)\n";
print Rinput "print(ks_test[1])\n";
print Rinput "print(ks_test[2])\n";
print Rinput "if (ks_test[2] >= 0.05){myKS = 'local dynamics of disease state does NOT differ significantly from ortholog'}\n";
print Rinput "if (ks_test[2] < 0.05){myKS = 'local dynamics of disease state IS significantly different from ortholog'}\n";
print Rinput "print(myKS)\n";
print Rinput "sink()\n";#quit
############
# barplot
print Rinput "d1 = data.frame(location=$location1, effect=factor($effect1), AAquery=factor($queryAA1),label=factor($label1), Yval=$KLdivergence1)\n";
print Rinput "d2 = data.frame(location=$location2, effect=factor($effect2), AAquery=factor($queryAA2),label=factor($label2), Yval=$KLdivergence2)\n";
print Rinput "myplot1 <- ggplot(data = d1, mapping = aes(x = location, y = Yval, fill=effect)) + xlim(0, $lengthID) + labs(x = 'orthologous mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity') + geom_text(aes(label = $label1), position = position_stack(vjust = 0.95), color = 'white', size = 1) + scale_fill_manual(values = c('darkgreen', 'darkred'))  + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey50'), panel.grid.major = element_line(colour = 'grey70'), panel.grid.minor = element_line(colour = 'grey70'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = location, y = Yval, fill=effect)) + xlim(0, $lengthID) + ggtitle(myKS) + labs(x = 'disease mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity') + geom_text(aes(label = $label2), position = position_stack(vjust = 0.95), color = 'white', size = 2) + scale_fill_manual(values = c('darkgreen', 'darkred')) + theme(axis.title.y = element_text(size=9), plot.title = element_text(color='blue', size=10, face='bold.italic'),legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey50'), panel.grid.major = element_line(colour = 'grey70'), panel.grid.minor = element_line(colour = 'grey70'))\n";
#print Rinput "myplot3 <- ggplot(data = d1, mapping = aes(x = location, y = Yval, fill=label)) + labs(x = 'orthologous mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
#print Rinput "myplot4 <- ggplot(data = d2, mapping = aes(x = location, y = Yval, fill=label)) + labs(x = 'disease mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(2, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
#print Rinput "print(myplot3, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
#print Rinput "print(myplot4, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
print Rinput "print(ks_test)\n";
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."KLmutations.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
print "THIS IS KS TEST FOR WHETHER MUTATION IN DISEASE AFFECTS MOLECULAR DYNAMICS SIGNIFICANTLY DIFFERENTLY THAN IN ORTHOLOG\n";
print "...result saved in file = KStest_diseaseVSortholog.txt\n\n";
print "close PDF viewer to continue\n";
print "done plotting bar chart of local impacts near points of mutation\n\n";
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.pdf\n";



# reread control files
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
      if ($header eq "homology"){$homology = $value;}
}
close IN;

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
my $seq_struct_flux = '';
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
    if ($header eq "test_type") { $seq_struct_flux = $value;}
	if ($header eq "color_scheme") { $color_scheme = $value;}

}
close IN;
sleep(1);

##########################################
if ($seq_struct_flux eq "flux"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomflux_ortholog/DROIDSfluctuation_$filenumber.txt") or next;
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
   print Rinput "data = read.table('atomflux_ortholog/DROIDSfluctuation_$filenumber.txt', header = TRUE)\n"; 
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
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
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
unlink("./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVGortholog.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGortholog.txt");
copy("myGranthamDistances_ortholog.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/myGranthamDistances.txt");
copy("mySeqSTATS_ortholog.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#copy("replylog.dat", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
}

#################### find average abs dFLUX #########################

open (OUT1, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGortholog.txt");
@dFLUXvals = ();
my @IN1 = <IN1>;
for (my $w = 0; $w < scalar @IN1; $w++){ # scan dFLUX
			     my $IN1row = $IN1[$w];
	             my @IN1row = split(/\s+/, $IN1row); 
			     my $dFLUXval = $IN1row[6]; # abs dFLUX
                 push (@dFLUXvals, $dFLUXval); 
                 }
    $statSCORE = new Statistics::Descriptive::Full; # avg abs delta flux
    $statSCORE->add_data (@dFLUXvals);
	$avg_abs_dFLUX = $statSCORE->mean();
    $avg_abs_dFLUX = sprintf "%.2f", $avg_abs_dFLUX;
close IN1;
print OUT1 "abs_dFLUX\t"."$avg_abs_dFLUX\n";
close OUT1;

#################### find overall RMSD #########################

#open (OUT2, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
#my @IN2 = <IN2>;
#for (my $w = 0; $w < scalar @IN2; $w++){ # scan for RMSD
#			     my $IN2row = $IN2[$w];
#	             my @IN2row = split(/\s+/, $IN2row); 
#			     $testheader = $IN2row[0];
#                 $RMSDtest = $IN2row[2];
#                 print OUT2 "test "."$testheader\t"."$RMSDtest\n";
#                 if ($testheader eq "Overall"){print OUT2 "overall_RMSD\t"."$RMSDtest\n";}
#                 }
#    
#close IN2;
#close OUT2;

################## adjust KStests.txt for multiple tests ############

sleep(1);
print " calculating multiple test correction \n\n";
sleep(1);
open (OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt") or die "could not create output file\n";
print OUT "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";
open (TMP, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt', header = TRUE)\n"; 
$posAA = "data\$posAA"; # position on reference structure	
$refAA = "data\$refAA"; # AA label on reference structure
$queryAA = "data\$queryAA"; # AA label on query structure
$Dval = "data\$Dval"; # D value for KS test
$pval = "data\$pval"; # p value for KS test
$signif = "data\$signif"; # significance label
print Rinput "p.adjust($pval, method = 'bonferroni', n = length($pval))\n";
print Rinput "p.adjust($pval, method = 'fdr', n = length($pval))\n"; #adjust p values for false discovery rate (i.e. Benjamini-Hochberg procedure)
print Rinput "write(p.adjust($pval, method = '$mtc', n = length($pval)), './DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt', sep='\n')\n";
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";	
close TMP;
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt") or die "could not create output file\n";
open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
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
print Rinput "library(gridExtra)\n";

# read data into R
if ($seq_struct_flux eq "flux"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGortholog.txt', header = TRUE)\n";} 
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt', header = TRUE)\n"; 
print Rinput "data3 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt', header = TRUE)\n"; 

if ($seq_struct_flux eq "flux"){
$pos_ref = "data1\$pos_ref"; # position on reference structure	
$res_ref = "data1\$res_ref"; # AA label on reference structure
$res_query = "data1\$res_query"; # AA label on query structure
$flux_ref_avg = "data1\$flux_ref_avg"; # avg flux on reference structure
$flux_query_avg = "data1\$flux_query_avg"; # avg flux on query structure
$delta_flux = "data1\$delta_flux"; # signed difference in avg flux
$abs_delta_flux = "data1\$abs_delta_flux"; # unsigned difference in avg flux
$delta_flux_kl = "data1\$KLdivergence"; # signed difference in flux as symmetric KL distance
}
$posAA = "data2\$posAA"; # position on reference structure	
$refAA = "data2\$refAA"; # AA label on reference structure
$queryAA = "data2\$queryAA"; # AA label on query structure
$Dval = "data2\$Dval"; # D value for KS test
$pval = "data2\$pval"; # p value for KS test
$signif = "data2\$signif"; # significance label

$label = "data3\$label"; # stat labels
$value = "data3\$value"; # stat values
print Rinput "data1\n";
print Rinput "data2\n";
print Rinput "data3\n";

# barplot
if ($seq_struct_flux eq "flux"){
print Rinput "d1A = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux_kl)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
print Rinput "mytable <- cbind(OVERALL=c('% AA match', 'avg Grantham distance', 'avg |dFLUX|'), STATISTICS=$value)\n";
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = '$refID')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = '$queryID')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'position (residue number)', y = 'dFLUX (signed KL distance)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dFLUX(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50')) +
  annotation_custom(tableGrob(mytable, theme = ttheme_minimal(base_size=8, base_colour='white')), xmin=0, ymin=0.2)\n";
#print Rinput "myplot4 <- ggplot() + annotation_custom(tableGrob(mytable, theme = ttheme_default(base_size=12), xmin=0, ymin=0))\n";
}

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
#print Rinput "print(myplot4, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSplot.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVGortholog.txt") or die "could not open file\n";
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
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSplot.pdf\n";
}

#############################################################################################

#################################################################

sub stats2 {

# run KS tests
print " running KS tests on each amino acid\n\n";
sleep(1);
print " reading control file\n\n";

my $queryID = '';
my $referenceID = '';
my $AA_count = '';
my $level_sig = '';
my $surface_or_ribbon = '';
my $seq_struct_flux = '';
my $color_scheme = '';


open(IN, "<"."DROIDSmut.ctl") or die "could not find CPPTRAJ input control file\n";
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
    if ($header eq "test_type") { $seq_struct_flux = $value;}
	if ($header eq "color_scheme") { $color_scheme = $value;}

}
close IN;
sleep(1);

##########################################
if ($seq_struct_flux eq "flux"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomflux_disease/DROIDSfluctuation_$filenumber.txt") or next;
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
   print Rinput "data = read.table('atomflux_disease/DROIDSfluctuation_$filenumber.txt', header = TRUE)\n"; 
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
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
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
unlink("./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVGdisease.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGdisease.txt");
copy("myGranthamDistances_disease.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/myGranthamDistances.txt");
copy("mySeqSTATS_disease.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#copy("replylog.dat", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
}

#################### find average abs dFLUX #########################

open (OUT1, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGdisease.txt");
@dFLUXvals = ();
my @IN1 = <IN1>;
for (my $w = 0; $w < scalar @IN1; $w++){ # scan dFLUX
			     my $IN1row = $IN1[$w];
	             my @IN1row = split(/\s+/, $IN1row); 
			     my $dFLUXval = $IN1row[6]; # abs dFLUX
                 push (@dFLUXvals, $dFLUXval); 
                 }
    $statSCORE = new Statistics::Descriptive::Full; # avg abs delta flux
    $statSCORE->add_data (@dFLUXvals);
	$avg_abs_dFLUX = $statSCORE->mean();
    $avg_abs_dFLUX = sprintf "%.2f", $avg_abs_dFLUX;
close IN1;
print OUT1 "abs_dFLUX\t"."$avg_abs_dFLUX\n";
close OUT1;

#################### find overall RMSD #########################

#open (OUT2, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
#my @IN2 = <IN2>;
#for (my $w = 0; $w < scalar @IN2; $w++){ # scan for RMSD
#			     my $IN2row = $IN2[$w];
#	             my @IN2row = split(/\s+/, $IN2row); 
#			     $testheader = $IN2row[0];
#                 $RMSDtest = $IN2row[2];
#                 print OUT2 "test "."$testheader\t"."$RMSDtest\n";
#                 if ($testheader eq "Overall"){print OUT2 "overall_RMSD\t"."$RMSDtest\n";}
#                 }
#    
#close IN2;
#close OUT2;

################## adjust KStests.txt for multiple tests ############

sleep(1);
print " calculating multiple test correction \n\n";
sleep(1);
open (OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt") or die "could not create output file\n";
print OUT "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";
open (TMP, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt', header = TRUE)\n"; 
$posAA = "data\$posAA"; # position on reference structure	
$refAA = "data\$refAA"; # AA label on reference structure
$queryAA = "data\$queryAA"; # AA label on query structure
$Dval = "data\$Dval"; # D value for KS test
$pval = "data\$pval"; # p value for KS test
$signif = "data\$signif"; # significance label
print Rinput "p.adjust($pval, method = 'bonferroni', n = length($pval))\n";
print Rinput "p.adjust($pval, method = 'fdr', n = length($pval))\n"; #adjust p values for false discovery rate (i.e. Benjamini-Hochberg procedure)
print Rinput "write(p.adjust($pval, method = '$mtc', n = length($pval)), './DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt', sep='\n')\n";
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";	
close TMP;
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt") or die "could not create output file\n";
open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
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
print Rinput "library(gridExtra)\n";

# read data into R
if ($seq_struct_flux eq "flux"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGdisease.txt', header = TRUE)\n";} 
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt', header = TRUE)\n"; 
print Rinput "data3 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt', header = TRUE)\n"; 

if ($seq_struct_flux eq "flux"){
$pos_ref = "data1\$pos_ref"; # position on reference structure	
$res_ref = "data1\$res_ref"; # AA label on reference structure
$res_query = "data1\$res_query"; # AA label on query structure
$flux_ref_avg = "data1\$flux_ref_avg"; # avg flux on reference structure
$flux_query_avg = "data1\$flux_query_avg"; # avg flux on query structure
$delta_flux = "data1\$delta_flux"; # signed difference in avg flux
$abs_delta_flux = "data1\$abs_delta_flux"; # unsigned difference in avg flux
$delta_flux_kl = "data1\$KLdivergence"; # signed difference in flux as symmetric KL distance
}
$posAA = "data2\$posAA"; # position on reference structure	
$refAA = "data2\$refAA"; # AA label on reference structure
$queryAA = "data2\$queryAA"; # AA label on query structure
$Dval = "data2\$Dval"; # D value for KS test
$pval = "data2\$pval"; # p value for KS test
$signif = "data2\$signif"; # significance label

$label = "data3\$label"; # stat labels
$value = "data3\$value"; # stat values
print Rinput "data1\n";
print Rinput "data2\n";
print Rinput "data3\n";

# barplot
if ($seq_struct_flux eq "flux"){
print Rinput "d1A = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux_kl)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
print Rinput "mytable <- cbind(OVERALL=c('% AA match', 'avg Grantham distance', 'avg |dFLUX|'), STATISTICS=$value)\n";
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = '$refID')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = '$queryID')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'position (residue number)', y = 'dFLUX (signed KL distance)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dFLUX(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50')) +
  annotation_custom(tableGrob(mytable, theme = ttheme_minimal(base_size=8, base_colour='white')), xmin=0, ymin=0.2)\n";
#print Rinput "myplot4 <- ggplot() + annotation_custom(tableGrob(mytable, theme = ttheme_default(base_size=12), xmin=0, ymin=0))\n";
}

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
#print Rinput "print(myplot4, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSplot.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVGdisease.txt") or die "could not open file\n";
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
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSplot.pdf\n";
}
#################################################################

sub stats3 {

# run KS tests
print " running KS tests on each amino acid\n\n";
sleep(1);
print " reading control file\n\n";

my $queryID = '';
my $referenceID = '';
my $AA_count = '';
my $level_sig = '';
my $surface_or_ribbon = '';
my $seq_struct_flux = '';
my $color_scheme = '';


open(IN, "<"."DROIDSnovel.ctl") or die "could not find CPPTRAJ input control file\n";
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
    if ($header eq "test_type") { $seq_struct_flux = $value;}
	if ($header eq "color_scheme") { $color_scheme = $value;}

}
close IN;
sleep(1);

##########################################
if ($seq_struct_flux eq "flux"){

mkdir ("DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $startN + $r;
   
   # collect AA info
    open(INFO, "<"."atomflux_novel/DROIDSfluctuation_$filenumber.txt") or next;
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
   print Rinput "data = read.table('atomflux_novel/DROIDSfluctuation_$filenumber.txt', header = TRUE)\n"; 
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
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
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
unlink("./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVGnovel.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGnovel.txt");
copy("myGranthamDistances_novel.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/myGranthamDistances.txt");
copy("mySeqSTATS_novel.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#copy("replylog.dat", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
}

#################### find average abs dFLUX #########################

open (OUT1, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGnovel.txt");
@dFLUXvals = ();
my @IN1 = <IN1>;
for (my $w = 0; $w < scalar @IN1; $w++){ # scan dFLUX
			     my $IN1row = $IN1[$w];
	             my @IN1row = split(/\s+/, $IN1row); 
			     my $dFLUXval = $IN1row[6]; # abs dFLUX
                 push (@dFLUXvals, $dFLUXval); 
                 }
    $statSCORE = new Statistics::Descriptive::Full; # avg abs delta flux
    $statSCORE->add_data (@dFLUXvals);
	$avg_abs_dFLUX = $statSCORE->mean();
    $avg_abs_dFLUX = sprintf "%.2f", $avg_abs_dFLUX;
close IN1;
print OUT1 "abs_dFLUX\t"."$avg_abs_dFLUX\n";
close OUT1;

#################### find overall RMSD #########################

#open (OUT2, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
#my @IN2 = <IN2>;
#for (my $w = 0; $w < scalar @IN2; $w++){ # scan for RMSD
#			     my $IN2row = $IN2[$w];
#	             my @IN2row = split(/\s+/, $IN2row); 
#			     $testheader = $IN2row[0];
#                 $RMSDtest = $IN2row[2];
#                 print OUT2 "test "."$testheader\t"."$RMSDtest\n";
#                 if ($testheader eq "Overall"){print OUT2 "overall_RMSD\t"."$RMSDtest\n";}
#                 }
#    
#close IN2;
#close OUT2;

################## adjust KStests.txt for multiple tests ############

sleep(1);
print " calculating multiple test correction \n\n";
sleep(1);
open (OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt") or die "could not create output file\n";
print OUT "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";
open (TMP, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt', header = TRUE)\n"; 
$posAA = "data\$posAA"; # position on reference structure	
$refAA = "data\$refAA"; # AA label on reference structure
$queryAA = "data\$queryAA"; # AA label on query structure
$Dval = "data\$Dval"; # D value for KS test
$pval = "data\$pval"; # p value for KS test
$signif = "data\$signif"; # significance label
print Rinput "p.adjust($pval, method = 'bonferroni', n = length($pval))\n";
print Rinput "p.adjust($pval, method = 'fdr', n = length($pval))\n"; #adjust p values for false discovery rate (i.e. Benjamini-Hochberg procedure)
print Rinput "write(p.adjust($pval, method = '$mtc', n = length($pval)), './DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt', sep='\n')\n";
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";	
close TMP;
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStests.txt") or die "could not create output file\n";
open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustPvalues.txt") or die "could not create output file\n";
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
print Rinput "library(gridExtra)\n";

# read data into R
if ($seq_struct_flux eq "flux"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGnovel.txt', header = TRUE)\n";} 
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/adjustKStests.txt', header = TRUE)\n"; 
print Rinput "data3 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt', header = TRUE)\n"; 

if ($seq_struct_flux eq "flux"){
$pos_ref = "data1\$pos_ref"; # position on reference structure	
$res_ref = "data1\$res_ref"; # AA label on reference structure
$res_query = "data1\$res_query"; # AA label on query structure
$flux_ref_avg = "data1\$flux_ref_avg"; # avg flux on reference structure
$flux_query_avg = "data1\$flux_query_avg"; # avg flux on query structure
$delta_flux = "data1\$delta_flux"; # signed difference in avg flux
$abs_delta_flux = "data1\$abs_delta_flux"; # unsigned difference in avg flux
$delta_flux_kl = "data1\$KLdivergence"; # signed difference in flux as symmetric KL distance
}
$posAA = "data2\$posAA"; # position on reference structure	
$refAA = "data2\$refAA"; # AA label on reference structure
$queryAA = "data2\$queryAA"; # AA label on query structure
$Dval = "data2\$Dval"; # D value for KS test
$pval = "data2\$pval"; # p value for KS test
$signif = "data2\$signif"; # significance label

$label = "data3\$label"; # stat labels
$value = "data3\$value"; # stat values
print Rinput "data1\n";
print Rinput "data2\n";
print Rinput "data3\n";

# barplot
if ($seq_struct_flux eq "flux"){
print Rinput "d1A = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux_kl)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
print Rinput "mytable <- cbind(OVERALL=c('% AA match', 'avg Grantham distance', 'avg |dFLUX|'), STATISTICS=$value)\n";
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = '$refID')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = '$queryID')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'position (residue number)', y = 'dFLUX (signed KL distance)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dFLUX(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50')) +
  annotation_custom(tableGrob(mytable, theme = ttheme_minimal(base_size=8, base_colour='white')), xmin=0, ymin=0.2)\n";
#print Rinput "myplot4 <- ggplot() + annotation_custom(tableGrob(mytable, theme = ttheme_default(base_size=12), xmin=0, ymin=0))\n";
}

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
#print Rinput "print(myplot4, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSplot.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVGnovel.txt") or die "could not open file\n";
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
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSplot.pdf\n";
}

#############################################################################################
sub novelMD {
#################################     
     # create mutate_protein.cmd script
open(MUT, ">"."mutate_list_novel.txt");
print MUT "substitution\t"."position\n";
close MUT;
print "opening mutate_list_novel.txt using gedit\n\n";
print "type two tab separated columns under 'substitution' and 'position' then save and close\n\n";
print "for example\n\n";
print "substitution\t"."position\n";
print "ALA\t"."23\n";
print "TYR\t"."31\n";
print "PRO\t"."35\n";
print "LEU\t"."47\n";
print "ARG\t"."52\n";

system "gedit mutate_list_novel.txt\n";

# run mutate_protein.cmd script
print "\nmutant list file (mutate_list_novel.txt) was created\n";
sleep(2);
##################################     
# create mutate_protein.cmd script
open (LST, "<"."mutate_list_novel.txt") || die "could not find mutate_list.txt\n";
@LST = <LST>;
open(MUT, ">"."mutate_protein.cmd");
print MUT "open $refID".".pdb\n";
    for (my $l = 0; $l < scalar @LST; $l++){
        if ($l == 0){next;}
        $LSTrow = $LST[$l];
        @LSTrow = split(/\s+/, $LSTrow);
        $subsTYPE = $LSTrow[0];
        $subsPOS = $LSTrow[1];
        print MUT "swapaa $subsTYPE"." #0:$subsPOS\n";
        }
print MUT "write 0 $refID"."_novel.pdb\n";
close MUT;

# run mutate_protein.cmd script
system("$chimera_path"."chimera --nogui mutate_protein.cmd > mutate_protein.log\n");
print "\nnovel mutant protein PDB file was created\n";     
sleep(2);
#####################################
#system "pdb4amber -i $refID"."_novel.pdb -o $refID"."_novelREDUCED.pdb --dry --reduce \n";
#sleep(2);
#####################################
system "perl GUI_STATSMD_DROIDSds1.pl\n";
}
#############################################################################################
sub novelSVM {

my $queryID = '';
my $referenceID = '';
my $AA_count = '';
my $level_sig = '';
my $surface_or_ribbon = '';
my $seq_struct_flux = '';
my $color_scheme = '';

# reading DROIDS control file

open(IN, "<"."DROIDSnovel.ctl") or die "could not find DROIDS control file\n";
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
    if ($header eq "test_type") { $seq_struct_flux = $value;}
	if ($header eq "color_scheme") { $color_scheme = $value;}

}
close IN;
sleep(1);

## setting SVM parameters
system "perl GUI_SVM_DROIDSds1.pl\n";
open(IN, "<"."DROIDSsvm.ctl") or die "could not find SVM control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "kernel") { $kernel = $value;}
    if ($header eq "degree") { $degree = $value;}
    if ($header eq "gamma") { $gamma = $value;}
    }
close IN;
sleep(1);

## creating mutate_flux_novel.txt file
print "collecting local atom fluctuations near points of novel mutation\n";
sleep(1);
open(OUT1, ">"."mutate_flux_novel.txt") or die "could not open output file\n";
print OUT1 "mut_number\t"."subsLOCATION\t"."posAA\t"."refAA\t"."queryAA\t"."atomTYPE\t"."flux_ref\t"."flux_query\n";
open(MUT, "<"."mutate_list_novel.txt") or die "could not find list of mutations\n";
my @MUT = <MUT>;
for (my $m = 0; $m <= scalar @MUT; $m++){
  my $MUTrow = $MUT[$m];
  my @MUTrow = split(/\s+/, $MUTrow);
  my $subsTYPE = $MUTrow[0];
  if ($m => 1){$subsLOCATION = $MUTrow[1]}
  my $mut_number = $m;
  if ($subsLOCATION =~ m/\d/){print "mutation is near "."$subsLOCATION\n"};
  for (my $r = 0; $r <= $lengthID; $r++){
   
   $filenumber = $startN + $r;
   open(INFO, "<"."atomflux_novel/DROIDSfluctuation_$filenumber.txt") or next;
   my @INFO = <INFO>;
   for (my $rr = 0; $rr <= scalar @INFO; $rr++){
   my $INFOrow = $INFO[$rr];
	        my @INFOrow = split(/\s+/, $INFOrow); 
	        my $posAA = $INFOrow[1];
		   my $refAA = $INFOrow[2];
		   my $queryAA = $INFOrow[3];
             my $atomTYPE = $INFOrow[5];
             my $flux_ref = $INFOrow[6];
		   my $flux_query = $INFOrow[7];
     #print "$subsLOCATION\t"."$posAA\n";
     #if ($subLOCATION - $posAA <= $mut_targetsize){print "$mut_number\t"."$posAA\t"."$refAA\t"."$atomTYPE\t"."$flux_ref\t"."$$flux_query\n"}
     if ($subsLOCATION =~ m/\d/ && abs($subsLOCATION - $posAA) <=  $mut_targetsize){print OUT1 "$mut_number\t"."$subsLOCATION\t"."$posAA\t"."$refAA\t"."$queryAA\t"."$atomTYPE\t"."$flux_ref\t"."$flux_query\n"}
    }
    close INFO;   
  }
}
close MUT;
close OUT1;     
print "done collecting local atom fluctuations near points of novel mutation\n\n";
sleep(1);

##collect data for SVM     
print "collecting data to train the SVM (support vector machine)\n";
sleep(1);
open (SVM, ">"."SVMtrain_data.txt") || die "could not open SVMtrain_data.txt\n";
print SVM "N_atom\t"."CA_atom\t"."C_atom\t"."O_atom\t"."class\n";
open (IN1, "<"."mutate_flux_ortholog.txt") || die "could not open mutate_flux_ortholog.txt\n";
@IN1 = <IN1>;
@mut_pos = ();
@Natom = ();
@CAatom = ();
@Catom = ();
@Oatom = ();
for (my $s = 0; $s < scalar @IN1; $s++){
        if ($s == 0){next;}
        $IN1row = $IN1[$s];
        @IN1row = split(/\s+/, $IN1row);
        $mut_pos = $IN1row[1];
        $atom_type = $IN1row[5];
        $flux_R = $IN1row[6];
        $flux_Q = $IN1row[7];
        $dFLUX = ($flux_Q - $flux_R);
        if ($atom_type eq "N"){push (@Natom, $dFLUX); push (@mut_pos, $mut_pos);}
        if ($atom_type eq "CA"){push (@CAatom, $dFLUX);}
        if ($atom_type eq "C"){push (@Catom, $dFLUX);}
        if ($atom_type eq "O"){push (@Oatom, $dFLUX);}
        }
$len_Natom = scalar @Natom;
$len_CAatom = scalar @CAatom;
$len_Catom = scalar @Catom;
$len_Oatom = scalar @Oatom;
$len_min_ortho = min($len_Natom, $len_CAatom, $len_Catom, $len_Oatom);

for (my $ss = 0; $ss < $len_min_ortho; $ss++){
     print SVM "$Natom[$ss]\t"."$CAatom[$ss]\t"."$Catom[$ss]\t"."$Oatom[$ss]\t"."ortholog\n";     
     }
open (IN2, "<"."mutate_flux_disease.txt") || die "could not open mutate_flux_disease.txt\n";
@IN2 = <IN2>;
@mut_pos = ();
@Natom = ();
@CAatom = ();
@Catom = ();
@Oatom = ();
for (my $s = 0; $s < scalar @IN2; $s++){
        if ($s == 0){next;}
        $IN2row = $IN2[$s];
        @IN2row = split(/\s+/, $IN2row);
        $mut_pos = $IN2row[1];
        $atom_type = $IN2row[5];
        $flux_R = $IN2row[6];
        $flux_Q = $IN2row[7];
        $dFLUX = ($flux_Q - $flux_R);
        if ($atom_type eq "N"){push (@Natom, $dFLUX); push (@mut_pos, $mut_pos);}
        if ($atom_type eq "CA"){push (@CAatom, $dFLUX);}
        if ($atom_type eq "C"){push (@Catom, $dFLUX);}
        if ($atom_type eq "O"){push (@Oatom, $dFLUX);}
        }
$len_Natom = scalar @Natom;
$len_CAatom = scalar @CAatom;
$len_Catom = scalar @Catom;
$len_Oatom = scalar @Oatom;
$len_min_disease = min($len_Natom, $len_CAatom, $len_Catom, $len_Oatom);

for (my $ss = 0; $ss < $len_min_disease; $ss++){
     print SVM "$Natom[$ss]\t"."$CAatom[$ss]\t"."$Catom[$ss]\t"."$Oatom[$ss]\t"."disease\n";     
     }
close SVM;
close IN1;
close IN2;

print "collecting data to deploy the SVM (support vector machine)\n";
sleep(1);
open (SVM, ">"."SVMdeploy_data.txt") || die "could not open SVMtrain_data.txt\n";
print SVM "N_atom\t"."CA_atom\t"."C_atom\t"."O_atom\t"."class\n";
open (IN1, "<"."mutate_flux_novel.txt") || die "could not open mutate_flux_novel.txt\n";
@IN1 = <IN1>;
$svm_size = scalar @IN1 - 1;
@mut_pos = ();
@Natom = ();
@CAatom = ();
@Catom = ();
@Oatom = ();
for (my $s = 0; $s < scalar @IN1; $s++){
        if ($s == 0){next;}
        $IN1row = $IN1[$s];
        @IN1row = split(/\s+/, $IN1row);
        $mut_pos = $IN1row[1];
        $atom_type = $IN1row[5];
        $flux_R = $IN1row[6];
        $flux_Q = $IN1row[7];
        $dFLUX = ($flux_Q - $flux_R);
        if ($atom_type eq "N"){push (@Natom, $dFLUX); push (@mut_pos, $mut_pos);}
        if ($atom_type eq "CA"){push (@CAatom, $dFLUX);}
        if ($atom_type eq "C"){push (@Catom, $dFLUX);}
        if ($atom_type eq "O"){push (@Oatom, $dFLUX);}
        }
$len_Natom = scalar @Natom;
$len_CAatom = scalar @CAatom;
$len_Catom = scalar @Catom;
$len_Oatom = scalar @Oatom;
$len_min_novel = min($len_Natom, $len_CAatom, $len_Catom, $len_Oatom);

for (my $ss = 0; $ss < $len_min_novel; $ss++){
     print SVM "$Natom[$ss]\t"."$CAatom[$ss]\t"."$Catom[$ss]\t"."$Oatom[$ss]\t"."novel\n";     
     }
close SVM;
close IN2;

## R pipe for SVM
print "running SVM classifications on mutations in novel variant\n";
sleep(2);
open (SVM1, ">"."SVMtest.txt") || die "could not open SVMtest.txt\n";
close SVM1;
open (SVM2, ">"."SVMpreds.txt") || die "could not open SVMpreds.txt\n";
close SVM2;
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
# load plotting libraries
print Rinput "library('e1071')\n";
print Rinput "library('ggplot2')\n";
print Rinput "library(gridExtra)\n";
print Rinput "library('dplyr')\n";
# read in data
print Rinput "data1 = read.table('SVMtrain_data.txt', header = TRUE)\n";
print Rinput "data1\n";
$class = "data1\$class"; # classifier
$Natom = "data1\$N_atom"; # N ataom
$CAatom = "data1\$CA_atom"; # CA atom
$Catom = "data1\$C_atom"; # C atom
$Oatom = "data1\$O_atom"; # O atom
print Rinput "dataframe1 = data.frame(N=$Natom, CA=$CAatom, C=$Catom, O=$Oatom, class=factor($class))\n";
print Rinput "subframe1 = data.frame(N=$Natom, CA=$CAatom, C=$Catom, O=$Oatom)\n";
print Rinput "print (dataframe1)\n";
print Rinput "data2 = read.table('SVMdeploy_data.txt', header = TRUE)\n";
print Rinput "data2\n";
$class = "data2\$class"; # classifier
$Natom = "data2\$N_atom"; # N ataom
$CAatom = "data2\$CA_atom"; # CA atom
$Catom = "data2\$C_atom"; # C atom
$Oatom = "data2\$O_atom"; # O atom
print Rinput "dataframe2 = data.frame(N=$Natom, CA=$CAatom, C=$Catom, O=$Oatom, class=factor($class))\n";
print Rinput "subframe2 = data.frame(N=$Natom, CA=$CAatom, C=$Catom, O=$Oatom)\n";
print Rinput "print (dataframe2)\n";
# run SVM
print Rinput "rnd <- sample(1:30, 1)\n";
print Rinput "set.seed(rnd)\n";
print Rinput "mysample = sample_n(dataframe1, $svm_size, replace=TRUE)\n";
print Rinput "print(mysample)\n";
#print Rinput "sample_subframe = mysample(SL, SW, PL, PW)\n";
#print Rinput "print(sample_subframe)\n";
print Rinput "print('RUNNING SVM in R   ...MAY TAKE SEVERAL MINUTES')\n";
if ($kernel eq "linear"){
print Rinput "mytest <- svm(class~N+CA+C+O, mysample, probability = FALSE, type = 'C-classification', kernel = '$kernel')\n";
}
if ($kernel eq "polynomial"){
print Rinput "mytest <- svm(class~N+CA+C+O, mysample, probability =TRUE, type = 'C-classification', kernel = '$kernel', degree = $degree)\n";
}
if ($kernel eq "radial"){
print Rinput "mytest <- svm(class~N+CA+C+O, mysample, probability =TRUE, type = 'C-classification', kernel = '$kernel', gamma = $gamma)\n";
}
print Rinput "mypred <- predict(mytest, subframe2, probability = FALSE, decision.values = TRUE)\n";
print Rinput "myplot1 <-ggplot(mysample, aes(N, O, colour = class)) + geom_point(size = 0.2) + ggtitle('ortho/disease training data') + labs(x = 'dFLUX on N atom', y = 'dFLUX on O atom') + scale_color_manual(values = c('coral', 'goldenrod')) + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <-ggplot(mysample, aes(CA, C, colour = class)) + geom_point(size = 0.2) + ggtitle('ortho/disease training data') + labs(x = 'dFLUX on CA atom', y = 'dFLUX on C atom') + scale_color_manual(values = c('coral', 'goldenrod')) + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";                                 
print Rinput "myplot3 <-ggplot(dataframe1, aes(N, O, colour = class)) + geom_point(size = 0.2) + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot4 <-ggplot(dataframe1, aes(CA, C, colour = class)) + geom_point(size = 0.2) + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";                                 
print Rinput "mySVM <- summary(mytest)\n";
print Rinput "print(mySVM)\n";
print Rinput "myplot5 <-ggplot(subframe2, aes(N, O, colour = mypred)) + geom_point(size = 0.2) + ggtitle('novel variant testing data') + scale_color_manual(values = c('coral', 'goldenrod')) + labs(x = 'dFLUX on N atom', y = 'dFLUX on O atom')  + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";                                 
print Rinput "myplot6 <-ggplot(subframe2, aes(CA, C, colour = mypred)) + geom_point(size = 0.2) + ggtitle('novel variant testing data') + scale_color_manual(values = c('coral', 'goldenrod')) + labs(x = 'dFLUX on CA atom', y = 'dFLUX on C atom') + theme(panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";                                 
print Rinput "library('grid')\n";
print Rinput "pushViewport(viewport(layout = grid.layout(2, 2)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))\n";
print Rinput "print(myplot5, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "print(myplot6, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))\n";
#print Rinput "print(myplot5, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
#print Rinput "print(myplot6, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))\n";
print Rinput "sink('SVMtest.txt')\n";
print Rinput "print (mySVM)\n";
print Rinput "sink()\n";
print Rinput "sink('SVMpreds.txt')\n";
print Rinput "print(mypred)\n";
print Rinput "sink()\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

## R plot for SVM     
print "creating R plots for mutations in novel variant\n";
sleep(1);     
print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$referenceID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSsvm.pdf";
copy($oldfilename, $newfilename);
my $oldfilename = "SVMpreds.txt";
my $newfilename = "./DROIDS_results_$queryID"."_$referenceID"."_$seq_struct_flux"."_$level_sig"."/"."SVMpreds.txt";
copy($oldfilename, $newfilename);
my $oldfilename = "SVMtest.txt";
my $newfilename = "./DROIDS_results_$queryID"."_$referenceID"."_$seq_struct_flux"."_$level_sig"."/"."SVMtest.txt";
copy($oldfilename, $newfilename);
sleep(1);
print " plotting subroutine is complete\n\n";
print " close PDF viewer to continue\n\n";
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSsvm.pdf\n";

print "SVM is completed\n\n";
sleep(1);

########################################################
# report final performance plot for novel mutations
########################################################

print "analyzing local atom fluctuations near points of novel mutation\n";
sleep(1);
open(IN, "<"."DROIDSnovel.ctl") or die "could not find CPPTRAJ input control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    print "$header\t"."$value\n";
    if ($header eq "query") { $queryID = $value;}
    if ($header eq "cutoff_value") { $level_sig = $value;}
    if ($header eq "test_type") { $seq_struct_flux = $value;}
close IN;
}


open(OUT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt")or die "could not open statistics.txt\n";
print OUT2 "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\n";
#open(OUT3, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations_novelANDortholog.txt")or die "could not open statistics.txt\n";
#print OUT3 "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\t"."class\n";

@fluxR = ();
@fluxQ = ();
open(MUT2, "<"."mutate_flux_novel.txt") or die "could not find list of mutations\n";
my @MUT2 = <MUT2>;
for (my $s = 0; $s <= scalar @MUT2; $s++){
  if ($s == 0){next;}
  my $MUTrow = $MUT2[$s];
  my $nextMUTrow = $MUT2[$s+1];
  my @MUTrow = split(/\s+/, $MUTrow);
  my @nextMUTrow = split(/\s+/, $nextMUTrow);
  my $mut_number = $MUTrow[0];
  my $next_number = $nextMUTrow[0];
  my $mut_location = $MUTrow[1];
  my $test_location = $MUTrow[2];
  my $ref_AA = $MUTrow[3];
  my $mut_AA = $MUTrow[4];
  my $mut_fluxR = $MUTrow[6];
  my $mut_fluxQ = $MUTrow[7];
  push (@fluxR, $mut_fluxR);
  push (@fluxQ, $mut_fluxQ);
  #print "$mut_number\t"."$next_number\n";
  if ($mut_location == $test_location){$label = "$ref_AA"."->"."$mut_AA"."$mut_location";}
  if ($mut_number != $next_number){
     print scalar @fluxR; print "\n";
     print scalar @fluxQ; print "\n";
     
     # calculate avg dFLUX
     $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
     $statSCORE->add_data (@fluxR);
	$flux_ref_avg = $statSCORE->mean();
     #$flux_ref_n = $statSCORE->count();
     #print "flux_ref_n\t"."$flux_ref_n\n";
	$statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
     $statSCORE->add_data (@fluxQ);
     $flux_query_avg = $statSCORE->mean();
     #$flux_query_n = $statSCORE->count();
     #print "flux_query_n\t"."$flux_query_n\n";
     $delta_flux = ($flux_ref_avg - $flux_query_avg);
     # calculate JS divergence
     open (TMP1, ">"."flux_mut_temp.txt") or die "could not create temp file\n";
     print TMP1 "flux_ref\t"."flux_query\n";
     for (my $t = 0; $t <= scalar @fluxQ; $t++){print TMP1 "$fluxR[$t]\t"; print TMP1 "$fluxQ[$t]\n";}
     close TMP1;
     open (TMP2, ">"."flux_mut_KL.txt") or die "could not create temp file\n";
     close TMP2;
     open (Rinput, "| R --vanilla")||die "could not start R command line\n";
     print Rinput "library('FNN')\n";
     print Rinput "data = read.table('flux_mut_temp.txt', header = TRUE)\n"; 
     $flux_ref = "data\$flux_ref"; # flux on reference residue
     $flux_query = "data\$flux_query"; # flux on query residue
     print Rinput "d1 = data.frame(fluxR=$flux_ref, fluxQ=$flux_query)\n";
     #print Rinput "print(d1)\n";
     print Rinput "myKL<-KL.dist($flux_ref, $flux_query, k=10)\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink('flux_mut_KL.txt')\n";
     print Rinput "print(myKL[10])\n";
     print Rinput "sink()\n";
     # write to output file and quit R
     print Rinput "q()\n";# quit R 
     print Rinput "n\n";# save workspace image?
     close Rinput;
     open (TMP3, "<"."flux_mut_KL.txt") or die "could not create temp file\n";
     my @TMP3 = <TMP3>;
     for (my $tt = 0; $tt <= scalar @TMP3; $tt++){
     $TMP3row = $TMP3[$tt];
     @TMP3row = split (/\s+/, $TMP3row);
     $header = $TMP3row[0];
     $value = $TMP3row[1];
     #print "$header\t"."$value\n";
     if ($header eq "[1]"){$KL = $value;}
     }
     if ($delta_flux >= 0){$effect = "PROTEIN_STABILIZED";} # make KL value negative if dFLUX is negative
     if ($delta_flux < 0){$effect = "PROTEIN_DESTABILIZED";} # make KL value negative if dFLUX is negative
     print "my KL is "."$KL\n";
     close TMP3;
     print OUT2 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\n";
     #print OUT3 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\t"."novel\n";
     @fluxR = ();
     @fluxQ = ();
     }
 
}
close MUT2;
close OUT2;
#close OUT3;
sleep(1);
print "\ndone analyzing local atom fluctuations near points of novel mutation\n\n";
sleep(1);

###############################
## collect SVM classifications
@SVMclass = ();
open(IN, "<"."SVMpreds.txt") or die "could not find SVMpreds.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $search = @INrow[0];
	 #print "search\t"."$search\n";
      if ($search eq "ortholog" || $search eq "disease" ){
          for (my $ii = 0; $ii < scalar @INrow; $ii++){
              $SVMclass = $INrow[$ii];
              #print "SVMclass\t"."$SVMclass\t";
              push(@SVMclass, $SVMclass);
              }
          }
      
}
close IN;
#print @SVMclass;
print "\n";
$size_svm_classarray = scalar @SVMclass;
print "size of svm class array = "."$size_svm_classarray\n";

open(IN, "<"."MDn.ctl") or die "could not find MDn.ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "Number_Runs"){$runsID = $value;}
      
}
close IN;
if ($runsID > 0){$step_svm_classarray = $size_svm_classarray/(8*$runsID);}
else {print "error - there are no sampling runs in MDn.ctl file\n"; exit;}
print "step size thru class array = "."$step_svm_classarray\n";

## create KLmutations_svm.txt##
open(OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations_svm.txt") or die "could not make output file\n";
print OUT "location\t"."startrefAA\t"."startqueryAA\t"."KLdivergence\t"."effect\t"."label\t"."percentDisease\t"."labelSVM\n";
open(IN, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt") or die "could not find input file\n";
my @IN = <IN>;
$location_svm_classarray = 0;
print "loc_mut\t"."location_svm_classarray\t"."percent_disease\t"."labelSVM\n";
for (my $i = 0; $i < scalar @IN; $i++){
	 if ($i == 0){next;}
      my @SVMortho = ();
      my @SVMdisease = ();
      my $INrow = $IN[$i];
      my @INrow = split (/\s+/, $INrow);
	 my $loc_mut = @INrow[0];
      chomp $INrow;
      $location_svm_classarray = int($loc_mut/$lengthID*$size_svm_classarray);
      $classSVM = $SVMclass[$location_svm_classarray];
      # collect classifications from nearby AA's to the mutation
      $collect_range = 2*$mut_targetsize + 1;
      for (my $ii = 0; $ii <= $collect_range; $ii++){
           $pointer = $location_svm_classarray - $mut_targetsize + $ii;
           $classSVM = $SVMclass[$pointer];
           if ($classSVM eq "ortholog"){push(@SVMortho, $classSVM);}
           if ($classSVM eq "disease"){push(@SVMdisease, $classSVM);}
           }
      #print @SVMortho;
      #print "\n";
      #print @SVMdisease;
      #print "\n";
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@SVMortho);
	 $SVMortho_count = $statSCORE->count();
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@SVMdisease);
	 $SVMdisease_count = $statSCORE->count();
      #print "$SVMdisease_count\n";
      #print "$SVMortho_count\n";
      $percent_disease = ($SVMdisease_count/($SVMdisease_count+$SVMortho_count))+0.001;
      if ($percent_disease <= 0.25){$labelSVM = "neutral";}
      if ($percent_disease > 0.25 && $percent_disease <= 0.5){$labelSVM = "likely_neutral";}
      if ($percent_disease > 0.5 && $percent_disease <= 0.75){$labelSVM = "likely_deleterious";}
      if ($percent_disease > 0.75){$labelSVM = "deleterious";}
      print "$loc_mut\t"."$location_svm_classarray\t"."$label\t"."$percent_disease\t"."$labelSVM\n";
      print OUT "$INrow\t"."$percent_disease\t"."$labelSVM\n";
      }
close IN;
close OUT;
sleep(1);

####################
print "\nplotting bar chart of local impacts near points of mutation\n\n";
sleep(1);
open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";

# read data into R
print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt', header = TRUE)\n"; 

$location1 = "data1\$location"; # position 	
$label1 = "data1\$label"; # AA label 
$effect1 = "data1\$effect"; # effect of mutation
$KLdivergence1 = "data1\$KLdivergence"; # impact of mutation
$queryAA1 = "data1\$startqueryAA"; # avg flux on query structure
print Rinput "data1\n";
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations_svm.txt', header = TRUE)\n"; 
$location2 = "data2\$location"; # position 	
$label2 = "data2\$label"; # AA label 
$effect2 = "data2\$effect"; # effect of mutation
$percentDisease = "data2\$percentDisease"; # percent of local region classified as disease dynamics
$labelSVM = "data2\$labelSVM"; # classification of novel mutation
$KLdivergence2 = "data2\$KLdivergence"; # impact of mutation
$queryAA2 = "data2\$startqueryAA"; # avg flux on query structure
print Rinput "data2\n";
########## global KS test #### does disease mutation significantly differ from orthologs
print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
print Rinput "print(ks_test)\n";
#print to file
#open(STAT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStest_novelVSortholog.txt")or die "could not open statistics.txt\n";
#close STAT;
#print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KStest_novelVSortholog.txt')\n";
#print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
#print Rinput "print(ks_test)\n";
#print Rinput "print(ks_test[1])\n";
#print Rinput "print(ks_test[2])\n";
#print Rinput "if (ks_test[2] >= 0.05){myKS = 'local dynamics of disease state does NOT differ significantly from ortholog'}\n";
#print Rinput "if (ks_test[2] < 0.05){myKS = 'local dynamics of disease state IS significantly different from ortholog'}\n";
#print Rinput "print(myKS)\n";
#print Rinput "sink()\n";#quit
print Rinput "myKS1 = 'SVM (support vector machine) classification of mutations in novel variant'\n";
print Rinput "myKS2 = 'percent of local flanking region around mutation with disease dynamics'\n";
############
# barplot
print Rinput "d1 = data.frame(location=$location1, effect=factor($effect1), AAquery=factor($queryAA1),label=factor($label1), Yval=$KLdivergence1)\n";
print Rinput "d2 = data.frame(location=$location2, labelSVM=factor($labelSVM), AAquery=factor($queryAA2),percentD = $percentDisease, label=factor($label2), Yval=$KLdivergence2)\n";
print Rinput "myplot1 <- ggplot(data = d1, mapping = aes(x = location, y = Yval, fill=effect)) + xlim(0, $lengthID) + labs(x = 'novel mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity') + geom_text(aes(label = $label1), position = position_stack(vjust = 0.95), color = 'white', size = 2) + scale_fill_manual(values = c('darkgreen', 'darkred')) + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey50'), panel.grid.major = element_line(colour = 'grey70'), panel.grid.minor = element_line(colour = 'grey70'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = location, y = Yval, fill=labelSVM)) + xlim(0, $lengthID) + ggtitle(myKS1) + labs(x = 'novel mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity') + geom_text(aes(label = $label2), position = position_stack(vjust = 0.95), color = 'gray75', size = 2) + scale_fill_brewer(palette='Accent') + theme(axis.title.y = element_text(size=9), plot.title = element_text(color='blue', size=10, face='bold.italic'),legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot3 <- ggplot(data = d2, mapping = aes(x = location, y = Yval, fill=percentD)) + xlim(0, $lengthID) + ggtitle(myKS2) + labs(x = 'novel mutation positions (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity') + geom_text(aes(label = $label2), position = position_stack(vjust = 0.95), color = 'gray75', size = 2) + scale_fill_gradient2(low = 'white', high = 'red', midpoint = 0.25) + theme(axis.title.y = element_text(size=9), plot.title = element_text(color='blue', size=10, face='bold.italic'),legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
print Rinput "print(myplot3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";

print Rinput "ks_test<-ks.test($KLdivergence1, $KLdivergence2)\n";
print Rinput "print(ks_test)\n";
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

print " copying plots\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."KLmutationsSVM.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
my $oldfilename = "DROIDSsvm.txt";
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSsvm.txt";
copy($oldfilename, $newfilename);	
sleep(1);
#print "THIS IS KS TEST FOR WHETHER MUTATION IN DISEASE AFFECTS MOLECULAR DYNAMICS SIGNIFICANTLY DIFFERENTLY THAN IN ORTHOLOG\n";
#print "...result saved in file = KStest_diseaseVSortholog.txt\n\n";
print "close PDF viewer to continue\n";
print "done plotting bar chart of local impacts near points of mutation\n\n";
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutationsSVM.pdf\n";



# reset control files
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
      if ($header eq "homology"){$homology = $value;}
}
close IN;


################################
# reread control files
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      
}
close IN;

print "final plotting is completed\n\n";
sleep(1);
}
#############################################################################################



