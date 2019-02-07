#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";
use Descriptive();
# specify the path to working directory for Chimera here

my $chimera_path = "/opt/UCSF/Chimera64-1.11/bin/";


#### Introductory Message #############

print "\n\nWelcome to Ambertools16 CPPTRAJ - for analyzing atom vector trajectories\n\n";
print "...and to DROIDS- Detecting Relative Outlier Impacts in Dynamic Simulations

- statistic procedures for functional evolutionary comparison
  of molecular dynamic simulation \n\n";

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### Declare variables ####
my @rep;
my $fileIDq = '';
my $fileIDr = '';
my $runsID = '';
my $implicit=0;
my $explicit=0;

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("CPPTRAJ atom vector trajectory analysis"); # Titles the main window

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
		
	my $runsFrame = $pdbFrame->Frame();
		my $runsLabel = $runsFrame->Label(-text=>"number of repeated MD runs: ");
		my $runsEntry = $runsFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$runsID
					);

# Solvation Frame
my $solnFrame = $mw->Frame(	-label => "Method of Solvation",
				-relief => "groove",
				-borderwidth => 2
				);
	my $implicitCheck = $solnFrame->Radiobutton( -text => "implicit",
						-value=>"im",
						-variable=>\$solvType
						);
	my $explicitCheck = $solnFrame->Radiobutton( -text => "explicit",
						-value=>"ex",
						-variable=>\$solvType
						);		

# Buttons
my $controlButton = $mw -> Button(-text => "create CPPTRJ .ctl files (first)", 
				-command => \&control
				); # Creates a file button

my $infoButton = $mw -> Button(-text => "create atom info files", 
				-command => \&info
				); # Creates a file button

my $fluxButton = $mw -> Button(-text => "create atom fluctuation files", 
				-command => \&flux
				); # Creates a file button

my $doneButton = $mw -> Button(-text => "parse / prepare files for DROIDS", 
				-command => \&done
				); # Creates a file button

my $stopButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop
				); # Creates a file button

#### Organize GUI Layout ####

$QfileLabel->pack(-side=>"left");
$QfileEntry->pack(-side=>"left");
$RfileLabel->pack(-side=>"left");
$RfileEntry->pack(-side=>"left");
$runsLabel->pack(-side=>"left");
$runsEntry->pack(-side=>"left");
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$doneButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$fluxButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$infoButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$controlButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
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

MainLoop; # Allows Window to Pop Up


########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################


sub stop {exit;}

########################################################################################

sub control { # Write a control file and then call appropriate scripts that reference control file
	
if ($solvType eq "im") {$implicit = 1;}
if ($solvType eq "ex") {$explicit = 1;}

### make atom info control files ###	
open(ctlFile1, '>', "atominfo_$fileIDq"."_0.ctl") or die "Could not open output file";
my $parm_label1 = '';
if ($implicit == 1) {my $parm_label1 = "vac_"."$fileIDq".".prmtop"; print ctlFile1 "parm $parm_label1\n";}
if ($explicit == 1) {my $parm_label1 = "wat_"."$fileIDq".".prmtop"; print ctlFile1 "parm $parm_label1\n";}
my $traj_label1 = "prod_"."$fileIDq"."_0".".nc";
print ctlFile1 "trajin $traj_label1\n";
print ctlFile1 "resinfo !(:WAT)\n"; # all residues but not water
print ctlFile1 "atominfo \@CA,C,O,N,H&!(:WAT)\n"; # mask for all protein backbone atoms eliminating water
close ctlFile1;

open(ctlFile2, '>', "atominfo_$fileIDr"."_0.ctl") or die "Could not open output file";
my $parm_label2 = '';
if ($implicit == 1) {my $parm_label2 = "vac_"."$fileIDr".".prmtop"; print ctlFile2 "parm $parm_label2\n";}
if ($explicit == 1) {my $parm_label2 = "wat_"."$fileIDr".".prmtop"; print ctlFile2 "parm $parm_label2\n";}
my $traj_label2 = "prod_"."$fileIDr"."_0".".nc";
print ctlFile2 "trajin $traj_label2\n";
print ctlFile2 "resinfo !(:WAT)\n"; # all residues but not water
print ctlFile2 "atominfo \@CA,C,O,N,H&!(:WAT)\n"; # mask for all protein backbone atoms eliminating water
close ctlFile2;



for (my $i = 0; $i < $runsID; $i++){
### make atom flux control files ###	
open(ctlFile3, '>', "atomflux_$fileIDq"."_$i.ctl") or die "Could not open output file";
my $parm_label3 = '';
if ($implicit == 1) {my $parm_label3 = "vac_"."$fileIDq".".prmtop"; print ctlFile3 "parm $parm_label3\n";}
if ($explicit == 1) {my $parm_label3 = "wat_"."$fileIDq".".prmtop"; print ctlFile3 "parm $parm_label3\n";}
my $traj_label3 = "prod_"."$fileIDq"."_$i".".nc";
print ctlFile3 "trajin $traj_label3\n";	
print ctlFile3 "rms first\n";
print ctlFile3 "average crdset MyAvg\n";
print ctlFile3 "run\n";
print ctlFile3 "rms ref MyAvg\n";
print ctlFile3 "atomicfluct out fluct_$fileIDq"."_$i.txt \@CA,C,O,N,H&!(:WAT)\n";
#print ctlFile3 "byatom\n"; # hash out for avg atom flux, unhash for total atom flux
print ctlFile3 "run\n";
close ctlFile3;

open(ctlFile4, '>', "atomflux_$fileIDr"."_$i.ctl") or die "Could not open output file";
my $parm_label4 = '';
if ($implicit == 1) {my $parm_label4 = "vac_"."$fileIDr".".prmtop"; print ctlFile4 "parm $parm_label4\n";}
if ($explicit == 1) {my $parm_label4 = "wat_"."$fileIDr".".prmtop"; print ctlFile4 "parm $parm_label4\n";}
my $traj_label4 = "prod_"."$fileIDr"."_$i".".nc";
print ctlFile4 "trajin $traj_label4\n";	
print ctlFile4 "rms first\n";
print ctlFile4 "average crdset MyAvg\n";
print ctlFile4 "run\n";
print ctlFile4 "rms ref MyAvg\n";
print ctlFile4 "atomicfluct out fluct_$fileIDr"."_$i.txt \@CA,C,O,N,H&!(:WAT)\n";
#print ctlFile4 "byatom\n";  # hash out for avg atom flux, unhash for total atom flux
print ctlFile4 "run\n";
close ctlFile4;

### make atom corr control files ###	
open(ctlFile5, '>', "atomcorr_$fileIDq"."_$i.ctl") or die "Could not open output file";
my $parm_label5 = '';
if ($implicit == 1) {my $parm_label5 = "vac_"."$fileIDq".".prmtop"; print ctlFile5 "parm $parm_label5\n";}
if ($explicit == 1) {my $parm_label5 = "wat_"."$fileIDq".".prmtop"; print ctlFile5 "parm $parm_label5\n";}
my $traj_label5 = "prod_"."$fileIDq"."_$i".".nc";
print ctlFile5 "trajin $traj_label5\n";	
#print ctlFile5 "atomiccorr out corr_$fileIDq"."_$i.txt \@CA,C,O,N,H&!(:WAT)\n";
print ctlFile5 "atomiccorr out corr_$fileIDq"."_$i.txt \!(:WAT)\n";
print ctlFile5 "run\n";
close ctlFile5;

open(ctlFile6, '>', "atomcorr_$fileIDr"."_$i.ctl") or die "Could not open output file";
my $parm_label6 = '';
if ($implicit == 1) {my $parm_label6 = "vac_"."$fileIDr".".prmtop"; print ctlFile6 "parm $parm_label6\n";}
if ($explicit == 1) {my $parm_label6 = "wat_"."$fileIDr".".prmtop"; print ctlFile6 "parm $parm_label6\n";}
my $traj_label6 = "prod_"."$fileIDr"."_$i".".nc";
print ctlFile6 "trajin $traj_label6\n";	
#print ctlFile6 "atomiccorr out corr_$fileIDr"."_$i.txt \@CA,C,O,N,H&!(:WAT)\n";
print ctlFile6 "atomiccorr out corr_$fileIDr"."_$i.txt \!(:WAT)\n";
print ctlFile6 "run\n";
close ctlFile6;

  } # end per run loop 
my $prefix = "";
open(metafile, '>', "$fileIDr.meta") or die "Could not open output file";
if ($implicit == 1) {$prefix = "vac";}
if ($explicit == 1) {$prefix = "wat";}
print metafile "amber\n$prefix"."_$fileIDr.prmtop\nprod_$fileIDr"."_0.nc\n";
close metafile;

print "\n\nall control files are done\n\n";


}

###################################################################################################

sub info { # launch atom info
system("cpptraj "."-i ./atominfo_$fileIDq"."_0.ctl | tee cpptraj_atominfo_$fileIDq.txt");
system("cpptraj "."-i ./atominfo_$fileIDr"."_0.ctl | tee cpptraj_atominfo_$fileIDr.txt");
}

###################################################################################################

sub flux { # launch atom fluctuation calc
for (my $i = 0; $i < $runsID; $i++){
system("cpptraj "."-i ./atomflux_$fileIDq"."_$i.ctl | tee cpptraj_atomflux_$fileIDq.txt");
system("cpptraj "."-i ./atomflux_$fileIDr"."_$i.ctl | tee cpptraj_atomflux_$fileIDr.txt");
  }
}

###################################################################################################

sub done {

print "Enter residue number at the start of both chains\n";
print "(e.g. enter 389 if starts at THR 389.A) \n";
print "(e.g. enter 1 if starts at MET 1.A) \n\n";
my $startN = <STDIN>;
chop($startN);

sleep(2);
print "\n\n searching for atom info file = "."cpptraj_atominfo_$fileIDr.txt\n";
sleep(2);
print "\n\n creating atom_residue_list_$fileIDr.txt\n";
open(OUT, ">"."atom_residue_list_$fileIDr.txt") or die "could open output file\n";
print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
open(IN, "<"."cpptraj_atominfo_$fileIDr.txt") or die "could not find atom info file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $atomnumber = @INrow[1];
	 my $atomlabel = @INrow[2];
	 my $resnumber = @INrow[3];
	 my $resindex = $resnumber + ($startN - 1);
	 my $reslabel = @INrow[4];
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resindex\t $reslabel\n"}
   }
close IN;
close OUT;
sleep(2);
print "\n\n creating atom_residue_list_$fileIDq.txt\n";
open(OUT, ">"."atom_residue_list_unmodified_$fileIDq.txt") or die "could open output file\n";
print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
open(IN, "<"."cpptraj_atominfo_$fileIDq.txt") or die "could not find atom info file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $atomnumber = @INrow[1];
	 my $atomlabel = @INrow[2];
	 my $resnumber = @INrow[3];
	 my $resindex = $resnumber + ($startN - 1);
	 my $reslabel = @INrow[4];
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resindex\t $reslabel\n"}
   }
close IN;
close OUT;
sleep(2);
#################################################################
print "\n\n modifying alignment for color mapping to reference structure\n";
sleep(2);  # need to remove indels from ref sequence and any corresponding AA's in query

open(OUT1, ">"."$fileIDr"."_alignREFMAP.aln") or die "could not open output file\n";
open(IN1, "<"."$fileIDr"."_align.aln") or die "could not open alignment...did you save as (.aln)?\n";
print OUT1 "CLUSTAL W ALN saved from UCSF Chimera/MultAlignViewer\n\n";
my @IN1 = <IN1>;
my $position = 0;
for (my $i = 0; $i < scalar @IN1; $i++){
	my $IN1row = $IN1[$i];
	my $IN1nextrow = $IN1[$i+1];
	if ($IN1row =~ m/$fileIDr/){my @IN1row = split(/\s+/, $IN1row); $header_ref = $IN1row[0]; $seq_ref =$IN1row[1]; print "$header_ref\t"."$seq_ref\n";
															my @IN1nextrow = split(/\s+/, $IN1nextrow); $header_query = $IN1nextrow[0]; $seq_query =$IN1nextrow[1]; print "$header_query\t"."$seq_query\n";
															my @seq_ref = split(//,$seq_ref);
															my @seq_query = split(//,$seq_query);
															my $new_seq_ref = "";
															my $new_seq_query = "";
															for (my $ii = 0; $ii < length $seq_ref; $ii++){
																      my $respos = $ii+1;
																			$position = $position+1;
																			my $AAref = @seq_ref[$ii]; 
																			my $AAquery = @seq_query[$ii];
																			if ($AAref ne "."){$new_seq_ref = $new_seq_ref.$AAref; $new_seq_query = $new_seq_query.$AAquery;}
															}
															print OUT1 "$header_ref\t"."$new_seq_ref\n";
															print OUT1 "$header_query\t"."$new_seq_query\n\n";
																													
															}
}
close OUT1;
close IN1;
sleep (2);

#################################################################
print "\n\n calculating AA sequence similarity and Grantham distances\n";
sleep(2);
print "\n\n creating vertical alignment files\n";
sleep(2);
open(OUT1, ">"."$fileIDr"."_vertalign_ref.aln") or die "could not open output file\n";
print OUT1 "respos\t"."seq_ref\n";
open(OUT2, ">"."$fileIDq"."_vertalign_query.aln") or die "could not open output file\n";
print OUT2 "respos\t"."seq_query\n";
open(OUT3, ">"."$fileIDr"."_vertalign_ref_indexed.aln") or die "could not open output file\n";
print OUT3 "respos\t"."seq_ref\n";
open(OUT4, ">"."$fileIDq"."_vertalign_query_indexed.aln") or die "could not open output file\n";
print OUT4 "respos\t"."seq_query\n";
open(OUT5, ">"."myGranthamDistances.txt") or die "could not open output file\n";
print OUT5 "respos\t"."distance\n";
open(IN1, "<"."$fileIDr"."_alignREFMAP.aln") or die "could not open alignment...did you save as (.aln)?\n";
my @IN1 = <IN1>;
my $position = 0;
my $positionINDEX = $startN-1;
my $AAsame_cnt = 0;
my $AA_cnt = 0;
my @gDISTS = ();
for (my $i = 0; $i < scalar @IN1; $i++){
	my $IN1row = $IN1[$i];
	my $IN1nextrow = $IN1[$i+1];
	if ($IN1row =~ m/$fileIDr/){my @IN1row = split(/\s+/, $IN1row); $header_ref = $IN1row[0]; $seq_ref =$IN1row[1]; print "$header_ref\t"."$seq_ref\n";
															my @IN1nextrow = split(/\s+/, $IN1nextrow); $header_query = $IN1nextrow[0]; $seq_query =$IN1nextrow[1]; print "$header_query\t"."$seq_query\n";
															my @seq_ref = split(//,$seq_ref);
															my @seq_query = split(//,$seq_query);
															for (my $ii = 0; $ii < length $seq_ref; $ii++){
																      $gDIST = '';
																			my $respos = $ii+1;
																			my $AAref = @seq_ref[$ii]; 
																			my $AAquery = @seq_query[$ii];
																			$position = $position+1;
																			$positionINDEX = $positionINDEX+1;
																			if ($AAquery eq $AAref){$AAsame_cnt = $AAsame_cnt + 1;}
																			if ($AAquery ne $AAref || $AAquery eq $AAref){$AA_cnt = $AA_cnt + 1;}
																			open(IN2, "<"."amino1to3.txt") or die "could not open amino1to3.txt\n";
																			my @IN2 = <IN2>;
																			#print "AAref "."$AAref\n";
																			#print "AAquery "."$AAquery\n";
																			for (my $iii = 0; $iii < scalar @IN2; $iii++){
																					my $AArow = @IN2[$iii];
																					my @AArow = split(/\s+/, $AArow);
																					$AAone = @AArow[0]; $AAthree = @AArow[1];
																					if ($AAone eq $AAref){print OUT1 "$position\t"."$AAthree\n"}
																					if ($AAone eq $AAquery){print OUT2 "$position\t"."$AAthree\n"}
																					if ($AAone eq $AAref){print OUT3 "$positionINDEX\t"."$AAthree\n"}
																					if ($AAone eq $AAquery){print OUT4 "$positionINDEX\t"."$AAthree\n"}
																			  	}
																			# determine Grantham distance
																			if ($AAquery eq $AAref || $AAquery eq "." || $AAref eq "."){$gDIST = 0;}
																			if ($AAquery ne $AAref && $AAquery ne "." && $AAref ne "."){
																			open(IN3, "<"."GranthamScores.txt") or die "could not open amino1to3.txt\n";
																			my @IN3 = <IN3>;
																			for (my $iiii = 0; $iiii < scalar @IN3; $iiii++){
																					my $GDrow = @IN3[$iiii];
																					my @GDrow = split(/\s+/, $GDrow);
																					$AAqueryTEST = @GDrow[0]; $AArefTEST = @GDrow[1]; $gDISTtest = @GDrow[2];
																					if(uc $AAqueryTEST eq $AAquery && uc $AArefTEST eq $AAref){$gDIST = $gDISTtest;} # grantham matrix
																					elsif(uc $AAqueryTEST eq $AAref && uc $AArefTEST eq $AAquery){$gDIST = $gDISTtest;} # to cover other half of matrix
																					}
																			}
																			
																			print OUT5 "$position\t"."$gDIST\n";
																			if ($gDIST > 0) {push (@gDISTS, $gDIST);} # average only non-zero Grantham Distances
																																																									 
													        }
															}
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close IN1;
close IN2;

### whole sequence stats
$statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@gDISTS);
					 $avg_gDIST = $statSCORE->mean();
					 $avg_gDIST = sprintf "%.2f", $avg_gDIST;

$AAseqsim = int(($AAsame_cnt/$AA_cnt+0.0001)*100);
$AAseq_matchfreq = ($AAsame_cnt/$AA_cnt+0.0001)*100;
$AAseq_matchfreq = sprintf "%.2f", $AAseq_matchfreq;
print "\n\nAA sequence similarity = "."$AAseqsim"."%\n";
print "avg Grantham Distance = "."$avg_gDIST"."\n";
open(OUT6, ">"."mySeqSTATS.txt") or die "could not open output file\n";
print OUT6 "label\t"."value\n";
print OUT6 "AAmatchFreq\t"."$AAseq_matchfreq\n";
print OUT6 "avgGranthamDist\t"."$avg_gDIST\n";
close OUT6;

sleep (2);


#################################################################################
print "\n\n adding gaps to query atom residue list (if needed)\n";
sleep(2);

open(IN1, "<"."atom_residue_list_unmodified_$fileIDq.txt") or die "could not open atom_residue_list.txt\n";
open(IN2, "<"."$fileIDq"."_vertalign_query_indexed.aln") or die "could not open output file\n";
open(OUT, ">"."atom_residue_list_$fileIDq.txt") or die "could not make atom_residue_list.txt\n";
#print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
my @IN1 = <IN1>;
my @IN2 = <IN2>;
$indelCount = 0;
for (my $i = 0; $i < scalar @IN1; $i++){ # scan residue list
	         my $IN1row = $IN1[$i];
			     my @IN1row = split(/\s+/, $IN1row);
			     my $atomnumber = $IN1row[0];
			     my $atomlabel = $IN1row[1];
			     my $resindex = $IN1row[2];
					 my $reslabel = $IN1row[3];
					 				 
					 for (my $j = 0; $j < scalar @IN2; $j++){ # scan alignment
			        my $IN2row = $IN2[$j];
			        my @IN2row = split(/\s+/, $IN2row);
			        my $pos_query = $IN2row[0] - $indelCount;
			        my $res_query = $IN2row[1];
							my $resindexgap = $resindex + $indelCount;
							#print "$pos_query\t"."$resindex\n";
				      if ($pos_query == $resindex && $res_query eq "xxx"){print OUT "na\t"."na\t"."na\t"."xxx\n"; print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n"; $indelCount = $indelCount+1}
							if ($pos_query == $resindex && $res_query ne "xxx"){print OUT "$atomnumber\t"."$atomlabel\t"."$resindexgap\t"."$reslabel\n";}		
			       
						 }
					 
					 #print "$reslabel\t".@skipped."\n";
					 
			}

close IN1;
close IN2;
close OUT;

#########################################################################################
##########  FLUX analysis     ###########################################################
#########################################################################################
print "\n\n collecting atomic fluctuation values (may take a minute)\n\n";
sleep(2);
open (OUT1, ">"."DROIDSfluctuation.txt") or die "could not create output file\n";
print OUT1 "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
open(IN3, "<"."atom_residue_list_$fileIDr.txt") or die "could not open atom_residue_list.txt\n";
open(IN4, "<"."atom_residue_list_$fileIDq.txt") or die "could not open atom_residue_list.txt\n";
my @IN3 = <IN3>;
my @IN4 = <IN4>;
      for (my $i = 0; $i < scalar @IN3; $i++){ # scan atom type
			     my $IN3row = $IN3[$i];
	         my @IN3row = split(/\s+/, $IN3row); 
			     my $atomnumberR = $IN3row[0]; my $atomlabelR = $IN3row[1]; my $resnumberR = $IN3row[2]; my $reslabelR = $IN3row[3];
					 my $IN4row = $IN4[$i];
	         my @IN4row = split(/\s+/, $IN4row); 
			     my $atomnumberQ = $IN4row[0]; my $atomlabelQ = $IN4row[1]; my $resnumberQ = $IN4row[2]; my $reslabelQ = $IN4row[3];
					 #print "atom+res REF"."$atomnumberR\t"."$atomlabelR\t"."$resnumberR\t"."$reslabelR\n";	                  
					 #print "atom+res QUERY"."$atomnumberQ\t"."$atomlabelQ\t"."$resnumberQ\t"."$reslabelQ\n";
					 # assemble fluctuation data
			     for (my $ii = 0; $ii < $runsID; $ii++){  #scan flux data
	            $sample = $ii;
							open(IN5, "<"."fluct_$fileIDq"."_$ii.txt") or die "could not open fluct file for $fileIDq\n";
              open(IN6, "<"."fluct_$fileIDr"."_$ii.txt") or die "could not open fluct file for $fileIDr\n";
	            my @IN5 = <IN5>;
              my @IN6 = <IN6>;
			        $flux_query = '';
							$flux_ref = '';
							for (my $iii = 0; $iii < scalar @IN5; $iii++){
							    my $IN5row = $IN5[$iii];
									$IN5row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN5row = split(/\s+/, $IN5row);
									my $Qtest_atom_decimal = $IN5row[0];
									my $Qtest_atom = int($Qtest_atom_decimal);
							    #print "Q "."$Qtest_atom\t"."$atomnumberQ\n";
									if($atomnumberQ eq $Qtest_atom){$flux_query = $IN5row[1];}
							  }	
							for (my $iii = 0; $iii < scalar @IN6; $iii++){
									my $IN6row = $IN6[$iii];
									$IN6row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN6row = split(/\s+/, $IN6row);
			            my $Rtest_atom_decimal = $IN6row[0];
									my $Rtest_atom = int($Rtest_atom_decimal);
									#print "R "."$Rtest_atom\t"."$atomnumberR\n";
									if($atomnumberR eq $Rtest_atom){$flux_ref = $IN6row[1];}
							  }
							
					    if($resnumberR =~/\d/ && $flux_query=~/\d/ && $reslabelQ ne "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    }
							if($resnumberR =~/\d/ && $reslabelQ eq "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."NA\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."NA\t"."NA\n";
							    }
							
							
							} 	
							
							close IN5;
              close IN6;
							
					
							}
					 
	
close IN3;
close IN4;
close OUT1;

########################################################################################
print "\n\n choose homology level for comparing backbone atom dynamics\n\n";
print " strict = collect only exact matching aligned residues\n";
print "          (e.g. position 5 -> LEU LEU)\n\n";
print " loose  = collect any aligned residues\n";
print "          (e.g. position 5 -> LEU LEU or position 5 -> LEU ALA)\n\n"; 
my $homology = <STDIN>;
chop($homology);

print "\n\n averaging DROIDSfluctuations by residue\n\n";
mkdir ("atomflux") or die "please delete atomflux folder from previous run\n";
open (IN, "<"."DROIDSfluctuation.txt") or die "could not create input file\n";
my @IN = <IN>;
open (OUT2, ">"."DROIDSfluctuationAVG.txt") or die "could not create output file\n";
print OUT2 "pos_ref\t"."res_ref\t"."res_query\t"."flux_ref_avg\t"."flux_query_avg\t"."delta_flux\t"."abs_delta_flux\n";
@REFfluxAvg = ();
@QUERYfluxAvg = ();
for (my $j = 0; $j < scalar @IN; $j++){ # scan atom type
			     my $INrow = $IN[$j];
	         my @INrow = split(/\s+/, $INrow); 
			     my $sample = $INrow[0];
					 my $pos_ref = $INrow[1];
					 my $res_ref = $INrow[2];
					 my $res_query = $INrow[3];
					 my $atomnumber = $INrow[4];
					 my $atomlabel = $INrow[5];
					 my $flux_ref = $INrow[6];
					 my $flux_query = $INrow[7];
					 push(@REFfluxAvg, $flux_ref);
					 push(@QUERYfluxAvg, $flux_query);
					 my $INnextrow = $IN[$j+1];
	         my @INnextrow = split(/\s+/, $INnextrow); 
			     my $next_pos = $INnextrow[1];
					 print OUT "$sample\t"."$pos_ref\t"."$res_ref\t"."$res_query\t"."$atomnumber\t"."$atomlabel\t"."$flux_ref\t"."$flux_query\n";
					 
					 if ($homology eq "loose"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_query ne "xxx"){  # loose homology = collect all aligned residues  
           open (OUT, ">"."./atomflux/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					 if ($pos_ref =~ m/\d/ && $j>1){
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFfluxAvg);
					 $flux_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYfluxAvg);
					 $flux_query_avg = $statSCORE->mean();
					 $delta_flux = ($flux_ref_avg - $flux_query_avg);
					 $abs_delta_flux = abs($flux_ref_avg - $flux_query_avg);
					 print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\n";
					 @REFfluxAvg = ();
           @QUERYfluxAvg = ();
					 }
					 if ($next_pos eq ''){next;}
					 }}
					 					 
					 if ($homology eq "strict"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_ref eq $res_query && $res_query ne "xxx"){ # strict homology = collect only exact matching residues  
           open (OUT, ">"."./atomflux/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					 if ($pos_ref =~ m/\d/ && $j>1){
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFfluxAvg);
					 $flux_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYfluxAvg);
					 $flux_query_avg = $statSCORE->mean();
					 $delta_flux = ($flux_ref_avg - $flux_query_avg);
					 $abs_delta_flux = abs($flux_ref_avg - $flux_query_avg);
					 print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\n";
					 @REFfluxAvg = ();
           @QUERYfluxAvg = ();
					 }
					 if ($next_pos eq ''){next;}
					 }}
					 
					 
																
}
close IN;
close OUT;
close OUT2;

sleep(2);

##################################################################################################
print "\n\n done parsing CPPTRAJ data files\n\n";
sleep(2);

system "perl GUI3_DROIDS.pl\n";	
}

##################################################################################################
