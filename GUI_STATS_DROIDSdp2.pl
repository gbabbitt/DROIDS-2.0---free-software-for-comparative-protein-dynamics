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
      #if ($header eq "start"){$startN = $value;}
      if ($header eq "homology"){$homology = $value;}
      if ($header eq "num_chains"){$chainN = $value;}
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

# chains Frame
my $chainFrame = $mw->Frame(	-label => "SELECT CHAIN TO ANALYZE",
				-relief => "groove",
				-borderwidth => 2
				);
	my $allRadio = $chainFrame->Radiobutton( -text => "all-for DROIDS image",
						-value=>"all",
						-variable=>\$chain
                              );
     $allRadio->pack();
################################						);
#  create list of chain labels
@alphabet = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
@chainlist = ();
if($chainN > 26) {print "warning - number of chains exceeds alphabet\n";}
for(my $l = 0; $l < $chainN; $l++){
     $letter = $alphabet[$l];
     push(@chainlist, $letter);
     }
print "chains in structure are...\n";
print @chainlist;
print "\n\n";
################################
     if($chainN > 1){
     for (my $i = 0; $i < scalar @chainlist; $i++){
     $chainRadio = "chainRadio"."$i";
     $chainRadio = $chainFrame->Radiobutton( -text => "chain $chainlist[$i]",
						-value=>"$chainlist[$i]",
						-variable=>\$chain
						);
     $chainRadio->pack();
     }
     }
		

# local mutation defn Frame				
my $mutFrame = $mw->Frame();
	my $MfileFrame = $mutFrame->Frame();
		my $MsizeLabel = $mutFrame->Label(-text=>"define local mutation effect size for rank analysis\n (e.g. 5 = analyze dFLUX within 5 AA's of subs site): ");
		my $MsizeEntry = $mutFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$mut_targetsize
					);
# local mutation defn Frame				
my $pathFrame = $mw->Frame();
	my $PfileFrame = $pathFrame->Frame();
		my $PLabel = $pathFrame->Label(-text=>"define path to previous MD results for normal binding\n (e.g. /home/greg/Desktop/DROIDS-2.0_1ytb_TF_DNA/): ");
		my $PEntry = $pathFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$normalpath
					);		
# Buttons

my $statsButton = $mw -> Button(-text => "make statistical comparisons and plots in R", 
				-command => \&stats,
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
my $rankButton = $mw -> Button(-text => "rank relative local impacts of each mutation", 
				-command => \&rank,
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
$chainFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$ctlButton->pack(-side=>"top",
			-anchor=>"s"
			);
$statsButton->pack(-side=>"top",
			-anchor=>"s"
			);
$mutFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$MsizeLabel->pack(-side=>"left");
$MsizeEntry->pack(-side=>"left");
$pathFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$PLabel->pack(-side=>"left");
$PEntry->pack(-side=>"left");
$rankButton->pack(-side=>"top",
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
sub next {system "perl GUI_IMAGE_DROIDSdp2.pl";}
########################################################################################
sub ctl {
# make control file for DROIDS	
print("Making ctl file...\n");

	$testStr = "flux"; $testStrLong = "fluctuation";  # file and folder labels
	
open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
print CTL "test_type\t"."$testStr\t # test method\n";
close CTL;
print("CTL file made\n");
}

########################################################################

sub rank{
#$mut_targetsize = 3;
print "collecting local atom fluctuations near points of mutation\n";
sleep(1);
open(OUT1, ">"."mutate_flux.txt") or die "could not open output file\n";
print OUT1 "mut_number\t"."subsLOCATION\t"."posAA\t"."refAA\t"."queryAA\t"."atomTYPE\t"."flux_ref\t"."flux_query\n";
open(MUT, "<"."mutate_list_trans_offset.txt") or die "could not find list of mutations\n";
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
   open(INFO, "<"."atomflux/DROIDSfluctuation_$filenumber.txt") or next;
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
print "done collecting local atom fluctuations near points of mutation\n\n";
sleep(1);

print "collecting local atom fluctuations from normally bound interaction\n";
sleep(1);
open(OUT1, ">"."mutate_flux_normal.txt") or die "could not open output file\n";
print OUT1 "mut_number\t"."subsLOCATION\t"."posAA\t"."refAA\t"."queryAA\t"."atomTYPE\t"."flux_ref\t"."flux_query\n";
open(MUT, "<"."mutate_list_trans_offset.txt") or die "could not find list of mutations\n";
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
   open(INFO, "<"."$normalpath"."atomflux/DROIDSfluctuation_$filenumber.txt") or next;
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
print "done collecting local atom fluctuations near points of mutation\n\n";
sleep(1);

print "analyzing local atom fluctuations near points of mutation\n";
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
print OUT2 "location\t"."refAA\t"."queryAA\t"."KLdivergence\t"."effect\t"."label\n";

@fluxR = ();
@fluxQ = ();
open(MUT2, "<"."mutate_flux.txt") or die "could not find list of mutations\n";
my @MUT2 = <MUT2>;
open(MUT3, "<"."mutate_flux_normal.txt") or die "could not find list of mutations\n"; # NOTE: normal bound protein is reference (not unbound protein)
my @MUT3 = <MUT3>;
for (my $s = 0; $s <= scalar @MUT2; $s++){
  if ($s == 0){next;}
  my $MUTrow = $MUT2[$s];
  my $nextMUTrow = $MUT2[$s+1];
  my $normMUTrow = $MUT3[$s]; # NOTE: normal bound protein is reference (not unbound protein)
  my @MUTrow = split(/\s+/, $MUTrow);
  my @nextMUTrow = split(/\s+/, $nextMUTrow);
  my @normMUTrow = split(/\s+/, $normMUTrow);
  my $mut_number = $MUTrow[0];
  my $next_number = $nextMUTrow[0];
  my $mut_location = $MUTrow[1];
  my $test_location = $MUTrow[2];
  my $ref_AA = $MUTrow[3];
  my $mut_AA = $MUTrow[4];
  my $mut_fluxR = $normMUTrow[6]; # NOTE: normal bound protein is reference (not unbound protein)
  my $mut_fluxQ = $MUTrow[6];
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
     if ($delta_flux >= 0){$effect = "BINDING_STABILIZED";} # make KL value negative if dFLUX is negative
     if ($delta_flux < 0){$effect = "BINDING_DESTABILIZED";} # make KL value negative if dFLUX is negative
     print "my KL is "."$KL\n";
     close TMP3;
     print OUT2 "$mut_location\t"."$ref_AA\t"."$mut_AA\t"."$KL\t"."$effect\t"."$label\n";
     @fluxR = ();
     @fluxQ = ();
     }
 
}
close MUT2;
close MUT3;
close OUT2;
sleep(1);
print "\ndone analyzing local atom fluctuations near points of mutation\n\n";
sleep(1);
print "\nplotting bar chart of local impacts near points of mutation\n\n";
sleep(1);
open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";

# read data into R
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.txt', header = TRUE)\n"; 

$location = "data\$location"; # position 	
$label = "data\$label"; # AA label 
$effect = "data\$effect"; # effect of mutation
$KLdivergence = "data\$KLdivergence"; # impact of mutation
$queryAA = "data\$queryAA"; # avg flux on query structure
print Rinput "data\n";

# barplot
print Rinput "d1 = data.frame(location=$location, effect=factor($effect), AAquery=factor($queryAA),label=factor($label), Yval=$KLdivergence)\n";
print Rinput "myplot1 <- ggplot(data = d1, mapping = aes(x = location, y = Yval, fill=effect)) + xlim(0, $lengthID+10) + labs(x = 'position (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <- ggplot(data = d1, mapping = aes(x = location, y = Yval, fill=label)) + xlim(0, $lengthID+10) + labs(x = 'position (residue number)', y = 'mutational impact = KL distance') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(4.0, 'mm'),legend.text = element_text(size = 8), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(2, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";

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
print " close PDF viewer to continue\n\n";
print "\ndone plotting bar chart of local impacts near points of mutation\n\n";
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/KLmutations.pdf\n";
}


#################################################################

sub stats {

print "Enter position of N terminal amino acid (default = 0)\n";
my $startN = <STDIN>;
chop($startN);
if ($startN eq ''){$startN = 0;}


if ($chain eq "all"){# stats subroutine for analyzing a single chain protein or all chains of muliti-chain protein
print " analyzing ONLY single chain protein or ALL chains as one\n\n";
sleep(2);

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
   
   $filenumber = $r;
   
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
copy("DROIDSfluctuationAVGchain.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGchain.txt");
copy("myGranthamDistances.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/myGranthamDistances.txt");
copy("mySeqSTATS.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
#copy("replylog.dat", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
}

#################### find average abs dFLUX #########################

open (OUT1, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/mySeqSTATS.txt");
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGchain.txt");
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
if ($seq_struct_flux eq "flux"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSfluctuationAVGchain.txt', header = TRUE)\n";} 
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
print Rinput "d1A = data.frame(pos1=$pos_ref+$startN, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref+$startN, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref+$startN, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux_kl)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref+$startN, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA+$startN, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
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
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/"."DROIDSplot_$chain.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVGchain.txt") or die "could not open file\n";
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
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/DROIDSplot_$chain.pdf\n";
}

#############################################################################################
#############################################################################################


if ($chain ne "all"){# stats subroutine for analyzing a specified chain (e.g. A, B or C)
print " analyzing chain $chain of multichain protein\n\n";
sleep(2);    

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

mkdir ("DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain") or die "DROID_results folder already exists...delete or rename it if running again";
open(STAT1, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStestsTEMP.txt")or die "could not open statistics.txt\n";
close STAT1;
open(STAT2, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStests.txt")or die "could not open statistics.txt\n";
print STAT2 "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\n";

for (my $r = 0; $r <= $AA_count; $r++){
   
   $filenumber = $r;
   
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
   print Rinput "sink('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStestsTEMP.txt')\n";
   print Rinput "ks_test<-ks.test($flux_ref, $flux_query)\n";
   print Rinput "print(ks_test)\n";
   print Rinput "sink()\n";#quit
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   print "\n\n";
   open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStestsTEMP.txt");
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
unlink("./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStestsTEMP.txt");
## copy residue avg data into results file
copy("DROIDSfluctuationAVGchain.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/DROIDSfluctuationAVGchain.txt");
copy("myGranthamDistances.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/myGranthamDistances.txt");
copy("mySeqSTATS.txt", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/mySeqSTATS.txt");
#copy("replylog.dat", "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."/replylog.dat");
}

#################### find average abs dFLUX #########################

open (OUT1, ">>"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/mySeqSTATS.txt");
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/DROIDSfluctuationAVGchain.txt");
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
open (OUT, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/adjustKStests.txt") or die "could not create output file\n";
print OUT "posAA\t"."refAA\t"."queryAA\t"."Dval\t"."pval\t"."signif\t"."chain\n";
open (TMP, ">"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/adjustPvalues.txt") or die "could not create output file\n";
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "data = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStests.txt', header = TRUE)\n"; 
$posAA = "data\$posAA"; # position on reference structure	
$refAA = "data\$refAA"; # AA label on reference structure
$queryAA = "data\$queryAA"; # AA label on query structure
$Dval = "data\$Dval"; # D value for KS test
$pval = "data\$pval"; # p value for KS test
$signif = "data\$signif"; # significance label
print Rinput "p.adjust($pval, method = 'bonferroni', n = length($pval))\n";
print Rinput "p.adjust($pval, method = 'fdr', n = length($pval))\n"; #adjust p values for false discovery rate (i.e. Benjamini-Hochberg procedure)
print Rinput "write(p.adjust($pval, method = '$mtc', n = length($pval)), './DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/adjustPvalues.txt', sep='\n')\n";
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";	
close TMP;
open (IN1, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/KStests.txt") or die "could not create output file\n";
open (IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/adjustPvalues.txt") or die "could not create output file\n";
open (IN3, "<"."./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/DROIDSfluctuationAVGchain.txt") or die "could not create output file\n";
my @IN1 = <IN1>;
my @IN2 = <IN2>;
my @IN3 = <IN3>;
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
                       # find chain ID
                       for (my $kkk = 0; $kkk < scalar @IN3; $kkk++){ # scan adjusted p values
				   my $IN3row = $IN3[$kkk]; my @IN3row = split(/\s+/, $IN3row);
				   my $posvalue = $IN3row[0];
                       my $chainvalue = $IN3row[8];
				   if ($posvalue == $posAA){$chain_ID = $chainvalue;}
				   }
                       
                   if ($posAA >= 1){print OUT "$posAA\t"."$refAA\t"."$queryAA\t"."$Dval\t"."$newPval\t"."$label_sig\t"."$chain_ID\n";}
				   
	
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
if ($seq_struct_flux eq "flux"){print Rinput "data1 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/DROIDSfluctuationAVGchain.txt', header = TRUE)\n";} 
print Rinput "data2 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/adjustKStests.txt', header = TRUE)\n"; 
print Rinput "data3 = read.table('./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/mySeqSTATS.txt', header = TRUE)\n"; 

# subset data for plotting a specific chain
print Rinput "data1chain <- subset(data1, chain=='$chain')\n"; 
#print Rinput "data2chain <- subset(data2, data2\$posAA==data1chain\$pos_ref)\n";
print Rinput "data2chain <- subset(data2, chain=='$chain')\n";

if ($seq_struct_flux eq "flux"){
$pos_ref = "data1chain\$pos_ref"; # position on reference structure	
$res_ref = "data1chain\$res_ref"; # AA label on reference structure
$res_query = "data1chain\$res_query"; # AA label on query structure
$flux_ref_avg = "data1chain\$flux_ref_avg"; # avg flux on reference structure
$flux_query_avg = "data1chain\$flux_query_avg"; # avg flux on query structure
$delta_flux = "data1chain\$delta_flux"; # signed difference in avg flux
$abs_delta_flux = "data1chain\$abs_delta_flux"; # unsigned difference in avg flux
$delta_flux_kl = "data1chain\$KLdivergence"; # signed difference in flux as symmetric KL distance

}
$posAA = "data2chain\$posAA"; # position on reference structure	
$refAA = "data2chain\$refAA"; # AA label on reference structure
$queryAA = "data2chain\$queryAA"; # AA label on query structure
$Dval = "data2chain\$Dval"; # D value for KS test
$pval = "data2chain\$pval"; # p value for KS test
$signif = "data2chain\$signif"; # significance label

$label = "data3\$label"; # stat labels
$value = "data3\$value"; # stat values

print Rinput "data1chain\n";
print Rinput "data2chain\n";
#print Rinput "data3\n";



# barplot
if ($seq_struct_flux eq "flux"){
print Rinput "d1A = data.frame(pos1=$pos_ref+$startN, label1r=$res_ref, label1q=$res_query, Y1val=$flux_ref_avg)\n";
print Rinput "d1B = data.frame(pos1=$pos_ref+$startN, label1r=$res_ref, label1q=$res_query, Y1val=$flux_query_avg)\n";
print Rinput "d2 = data.frame(pos2=$pos_ref+$startN, label2r=$res_ref, label2q=$res_query, Y2val=$delta_flux_kl)\n";
print Rinput "d3 = data.frame(pos3=$pos_ref+$startN, label3r=$res_ref, label3q=$res_query, Y3val=$abs_delta_flux)\n";
print Rinput "d4 = data.frame(pos4=$posAA+$startN, label4r=$refAA, label4q=$queryAA, Y4val = $Dval, Y4sig = $pval, Y4label = $signif)\n";;
#print Rinput "mytable <- cbind(OVERALL=c('% AA match', 'avg Grantham distance', 'avg |dFLUX|'), STATISTICS=$value)\n";
print Rinput "myplot1 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = d1A, mapping = aes(x = pos1, y = Y1val, color = '$refID')) + geom_line(data = d1B, mapping = aes(x = pos1, y = Y1val, color = '$queryID')) + theme(axis.title.y = element_text(size=9), axis.title.x=element_blank(), legend.title=element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
print Rinput "myplot2 <- ggplot(data = d2, mapping = aes(x = pos2, y = Y2val, fill=label2r)) + labs(x = 'chain $chain position (residue number)', y = 'dFLUX (signed KL distance)') + geom_bar(stat='identity')+ theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'),legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
#print Rinput "myplot3 <- ggplot(data = d3, mapping = aes(x = pos3, y = Y3val, fill=label3r)) + labs(x = 'position (residue number)', y = 'magnitude dFLUX(unsigned)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5), panel.background = element_rect(fill = 'grey30'))\n";
print Rinput "myplot3 <- ggplot(data = d4, mapping = aes(x = pos4, y = Y4val, fill=Y4label)) + labs(x = 'chain $chain position (residue number)', y = 'D value (KS test)') + geom_bar(stat='identity') + theme(axis.title.y = element_text(size=9), legend.title = element_blank(), panel.background = element_rect(fill = 'grey30'), panel.grid.major = element_line(colour = 'grey50'), panel.grid.minor = element_line(colour = 'grey50'))\n";
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
my $newfilename = "./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/"."DROIDSplot_$chain.pdf";
copy($oldfilename, $newfilename);	
sleep(1);
if ($scalingType eq "relative"){
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVGchain.txt") or die "could not open file\n";
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
system "evince ./DROIDS_results_$queryID"."_$refID"."_$seq_struct_flux"."_$level_sig"."_$chain"."/DROIDSplot_$chain.pdf\n";     
     
}

} # end sub

#############################################################################################