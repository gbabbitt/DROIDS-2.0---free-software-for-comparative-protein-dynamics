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

my $corrButton = $mw -> Button(-text => "create atom correlation files", 
				-command => \&corr
				); # Creates a file button

my $doneButton = $mw -> Button(-text => "align structures / prepare files for DROIDS", 
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
$corrButton->pack(-side=>"bottom",
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

sub corr { # launch atom correlation calc
for (my $i = 0; $i < $runsID; $i++){
system("cpptraj "."-i ./atomcorr_$fileIDq"."_$i.ctl | tee cpptraj_atomcorr_$fileIDq.txt");
system("cpptraj "."-i ./atomcorr_$fileIDr"."_$i.ctl | tee cpptraj_atomcorr_$fileIDr.txt");
  }
}

###################################################################################################

sub done {
	
print "IMPORTANT - Here you will need to run MatchMaker then Match-Align in Chimera\n\n";
print "            if satisfied with alignment, save as a clustal file with ref PDB ID\n";
print "            in title (e.g. 1ubq_align.aln)\n\n";
print "            DROIDS will prompt you to continue after closing Chimera\n\n";
print "continue? (y/n)\n";
my $go = <STDIN>;
chop($go);
if ($go eq "n") {exit;}
sleep(1);
print "            opening USCF Chimera and loading PDB ref structure\n\n";
print "            CREATE YOUR STRUCTURAL/SEQUENCE ALIGNMENT (.aln) NOW \n\n";
system("$chimera_path"."chimera $fileIDr.pdb $fileIDq.pdb\n");
sleep(2);
print "\n\n searching for atom info file = "."cpptraj_atominfo_$fileIDr.txt\n";
sleep(2);

###################################################################
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
	 my $reslabel = @INrow[4];
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resnumber\t $reslabel\n"}
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
	 my $reslabel = @INrow[4];
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resnumber\t $reslabel\n"}
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
print "\n\n creating vertical alignment files\n";
sleep(2);
open(OUT1, ">"."$fileIDr"."_vertalign_ref.aln") or die "could not open output file\n";
print OUT1 "respos\t"."seq_ref\n";
open(OUT2, ">"."$fileIDq"."_vertalign_query.aln") or die "could not open output file\n";
print OUT2 "respos\t"."seq_query\n";
open(IN1, "<"."$fileIDr"."_alignREFMAP.aln") or die "could not open alignment...did you save as (.aln)?\n";
my @IN1 = <IN1>;
my $position = 0;
for (my $i = 0; $i < scalar @IN1; $i++){
	my $IN1row = $IN1[$i];
	my $IN1nextrow = $IN1[$i+1];
	if ($IN1row =~ m/$fileIDr/){my @IN1row = split(/\s+/, $IN1row); $header_ref = $IN1row[0]; $seq_ref =$IN1row[1]; print "$header_ref\t"."$seq_ref\n";
															my @IN1nextrow = split(/\s+/, $IN1nextrow); $header_query = $IN1nextrow[0]; $seq_query =$IN1nextrow[1]; print "$header_query\t"."$seq_query\n";
															my @seq_ref = split(//,$seq_ref);
															my @seq_query = split(//,$seq_query);
															for (my $ii = 0; $ii < length $seq_ref; $ii++){
																      my $respos = $ii+1;
																			my $AAref = @seq_ref[$ii]; 
																			my $AAquery = @seq_query[$ii];
																			$position = $position+1;
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
																					
																			}
																																						 
													        }
															}
}
close OUT1;
close OUT2;
close IN1;
close IN2;
sleep (2);


#################################################################################
print "\n\n adding gaps to query atom residue list (if needed)\n";
sleep(2);

open(IN1, "<"."atom_residue_list_unmodified_$fileIDq.txt") or die "could not open atom_residue_list.txt\n";
open(IN2, "<"."$fileIDq"."_vertalign_query.aln") or die "could not open output file\n";
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
			     my $resnumber = $IN1row[2];
			     my $reslabel = $IN1row[3];
					 					 
					 for (my $j = 0; $j < scalar @IN2; $j++){ # scan alignment
			        my $IN2row = $IN2[$j];
			        my @IN2row = split(/\s+/, $IN2row);
			        my $pos_query = $IN2row[0] - $indelCount;
			        my $res_query = $IN2row[1];
				      if ($pos_query == $resnumber && $res_query eq "xxx"){print OUT "na\t"."na\t"."na\t"."xxx\n"; print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n"; $indelCount = $indelCount+1}
							if ($pos_query == $resnumber && $res_query ne "xxx"){print OUT "$atomnumber\t"."$atomlabel\t"."$resnumber\t"."$reslabel\n";}		
			       
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
							open(IN5, "<"."fluct_$fileIDr"."_$ii.txt") or die "could not open fluct file for $fileIDq\n";
              open(IN6, "<"."fluct_$fileIDq"."_$ii.txt") or die "could not open fluct file for $fileIDr\n";
	            my @IN5 = <IN5>;
              my @IN6 = <IN6>;
			        for (my $iii = 0; $iii < scalar @IN5; $iii++){
							    my $IN5row = $IN5[$iii];
									$IN5row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN5row = split(/\s+/, $IN5row);
									my $Qtest_atom_decimal = $IN5row[0];
									my $Qtest_atom = int($Qtest_atom_decimal);
							    #print "Q "."$Qtest_atom\t"."$atomnumberQ\n";
									if($atomnumberQ == $Qtest_atom){$flux_query = $IN5row[1];}
							  }	
							for (my $iii = 0; $iii < scalar @IN6; $iii++){
									my $IN6row = $IN6[$iii];
									$IN6row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN6row = split(/\s+/, $IN6row);
			            my $Rtest_atom_decimal = $IN6row[0];
									my $Rtest_atom = int($Rtest_atom_decimal);
									#print "R "."$Rtest_atom\t"."$atomnumberR\n";
									if($atomnumberR == $Rtest_atom){$flux_ref = $IN6row[1];}
							  }
							
					    if($resnumberR =~/\d/ && $flux_query=~/\d/ && $reslabelQ ne "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    }
							if($resnumberR =~/\d/ && $reslabelQ eq "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."NA\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."NA\n";
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
my @REFfluxAvg = ();
my @QUERYfluxAvg = ();
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
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFfluxAvg);
					 $flux_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYfluxAvg);
					 $flux_query_avg = $statSCORE->mean();
					 $delta_flux = -($flux_ref_avg - $flux_query_avg);
					 $abs_delta_flux = abs($flux_ref_avg - $flux_query_avg);
					 if ($pos_ref =~ m/\d/){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\n";}
					 my @REFfluxAvg = ();
           my @QUERYfluxAvg = ();
					 if ($next_pos eq ''){next;}
					 }}
					 					 
					 if ($homology eq "strict"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_ref eq $res_query && $res_query ne "xxx"){ # strict homology = collect only exact matching residues  
           open (OUT, ">"."./atomflux/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFfluxAvg);
					 $flux_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYfluxAvg);
					 $flux_query_avg = $statSCORE->mean();
					 $delta_flux = -($flux_ref_avg - $flux_query_avg);
					 $abs_delta_flux = abs($flux_ref_avg - $flux_query_avg);
					 if ($pos_ref =~ m/\d/){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\n";}
					 my @REFfluxAvg = ();
           my @QUERYfluxAvg = ();
					 if ($next_pos eq ''){next;}
					 }}
					 
					 
																
}
close IN;
close OUT;
close OUT2;

sleep(1);
#######################################################################################
print "  creating query FLUX scaled to avg reference FLUX\n";
sleep(2);
# find scaling factor
@vals = ();
open(IN, "<"."DROIDSfluctuationAVG.txt") or die "could not open file\n";
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
# created scaled data table
open(IN, "<"."DROIDSfluctuationAVG.txt") or die "could not open file\n";
open(OUT, ">"."DROIDSfluctuationAVGscaled.txt")or die "could not create scaled data\n";
print OUT "pos_ref\t"."res_ref\t"."res_query\t"."flux_ref_avg\t"."flux_query_avg\t"."delta_flux\t"."abs_delta_flux\n";
my @IN = <IN>;
for (my $c = 0; $c < scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $pos_ref = $INrow[0];
    my $res_ref = $INrow[1];
    my $res_query = $INrow[2];
    my $flux_ref_avg = $INrow[3];
    $flux_query_avg = $INrow[4] - $scaling_factor;
    $delta_flux = -($flux_ref_avg - $flux_query_avg);
    $abs_delta_flux = abs($delta_flux);
    if ($pos_ref ne "pos_ref"){print OUT "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\n";}
}
close IN;
close OUT;

print "\n\n parsing scaled DROIDSfluctuations by residue\n\n";
mkdir ("atomfluxscaled") or die "please delete atomflux folder from previous run\n";
open (IN, "<"."DROIDSfluctuation.txt") or die "could not create input file\n";
my @IN = <IN>;
my @REFfluxAvg = ();
my @QUERYfluxAvg = ();
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
					 $flux_query = $flux_query - $scaling_factor;
					 push(@REFfluxAvg, $flux_ref);
					 push(@QUERYfluxAvg, $flux_query);
					 my $INnextrow = $IN[$j+1];
	         my @INnextrow = split(/\s+/, $INnextrow); 
			     my $next_pos = $INnextrow[1];
					 print OUT "$sample\t"."$pos_ref\t"."$res_ref\t"."$res_query\t"."$atomnumber\t"."$atomlabel\t"."$flux_ref\t"."$flux_query\n";
					 if ($homology eq "loose"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_query ne "xxx"){  # loose homology = collect all aligned residues  
           open (OUT, ">"."./atomfluxscaled/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					 if ($next_pos eq ''){next;}
					 }}
					 if ($homology eq "strict"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_ref eq $res_query && $res_query ne "xxx"){ # strict homology = collect only exact matching residues  
           open (OUT, ">"."./atomfluxscaled/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					 if ($next_pos eq ''){next;}
					 }}
					 
					 
																
}
close IN;
close OUT;

#########################################################################################
##########  CORR analysis     ###########################################################
#########################################################################################
print "\n\n searching for atom info file = "."cpptraj_atominfo_$fileIDr.txt\n";
sleep(2);
open(IN, "<"."cpptraj_atominfo_$fileIDr.txt") or die "could not find atom info file\n";
my @IN = <IN>;
print @IN;
close IN;
print "\n\n CHOOSE residue name/number to reference all atom correlations
(e.g. 'leu15' = leucine at 15th residue or 'skip' to avoid CORR analyses)\n\n";
print " NOTE: this residue should typically coincide with a site of interest\n";
print "       such as site of mutation or site of important dynamic function\n\n";
sleep(1);
print " Chimera will now open the PBD of the reference structure so you can examine what amino acid to choose\n\n";
print " CLOSE CHIMERA TO CONTINUE\n\n";
system("$chimera_path"."chimera $fileIDr.pdb\n");
print " ENTER RESIDUE (e.g. leu15 or skip)\n";
my $res = <STDIN>;
chop($res);
sleep(2);
if ($res eq "skip"){goto SKIP;}
print "\n\n trimming atomic correlation files to $res\n\n";
$residue_number = substr($res, 3, length($res));
$amino = substr(uc $res, 0, 3);
print " labels are... "."$res\t"."$amino\t"."$residue_number\n\n";
print " this may take several minutes\n\n";

for (my $k = 0; $k < $runsID; $k++){  #scan corr data
	            print "trimming "."corr_$fileIDq"."_$k.txt"." and "."corr_$fileIDr"."_$k.txt\n"; 
							open(IN1, "<"."corr_$fileIDq"."_$k.txt") or die "could not open corr file for $fileIDq\n";
              open(OUT1, ">"."corr_trim_$fileIDq"."_$k.txt") or die "could not open corr output for $fileIDq\n";
							print OUT1 "resnumber\t"."atomnumber\t"."corr\n";
							open(IN2, "<"."corr_$fileIDr"."_$k.txt") or die "could not open corr file for $fileIDr\n";
              open(OUT2, ">"."corr_trim_$fileIDr"."_$k.txt") or die "could not open corr output for $fileIDr\n";
							print OUT2 "resnumber\t"."atomnumber\t"."corr\n";
							my @IN1 = <IN1>;
							my @IN2 = <IN2>;
							for (my $jj = 0; $jj < scalar @IN1; $jj++){ # scan corr files
								my $IN1row = $IN1[$jj];
								$IN1row =~ s/^\s+//;# need trim leading whitespace if present 
	              my @IN1row = split(/\s+/, $IN1row); 
			          my $pos1 = $IN1row[0];
								my $pos1 = int($pos1);
								my $test_resno1 = $IN1row[1];
								my $test_resno1 = int($test_resno1);
								my $corr1 = $IN1row[2];
								my $IN2row = $IN2[$jj];
								$IN2row =~ s/^\s+//;# need trim leading whitespace if present 
	              my @IN2row = split(/\s+/, $IN2row); 
			          my $pos2 = $IN2row[0];
								my $pos2 = int($pos2);
								my $test_resno2 = $IN2row[1];
								my $test_resno2 = int($test_resno2);
								my $corr2 = $IN2row[2];
								if ($residue_number == $test_resno1){print OUT1 "$test_resno1\t"."$pos1\t"."$corr1\n";}
								if ($residue_number == $test_resno2){print OUT2 "$test_resno2\t"."$pos2\t"."$corr2\n";}
								}
							  close IN1;
								close OUT1;
								close IN2;
								close OUT2;
			
       }


########################################################################################
			 
print "\n\n masking atomic correlation files to backbone atoms (CA, C, O, N)\n\n";			 
sleep(2);

for (my $k = 0; $k < $runsID; $k++){  #scan corr data
print "masking "."corr_trim_$fileIDq"."_$k.txt"." and "."corr_trim_$fileIDr"."_$k.txt\n";
open(OUT1, ">"."corr_mask_$fileIDq"."_$k.txt") or die "could not open corr output for $fileIDq\n";
open(OUT2, ">"."corr_mask_$fileIDr"."_$k.txt") or die "could not open corr output for $fileIDr\n";
open(IN1, "<"."corr_trim_$fileIDq"."_$k.txt") or die "could not open corr input for $fileIDq\n";
open(IN2, "<"."corr_trim_$fileIDr"."_$k.txt") or die "could not open corr input for $fileIDr\n";
open(IN3, "<"."atom_residue_list_$fileIDq.txt") or die "could not open atom residue list\n";
open(IN4, "<"."atom_residue_list_$fileIDr.txt") or die "could not open atom residue list\n";
my @IN1 = <IN1>;
my @IN2 = <IN2>;
my @IN3 = <IN3>;
my @IN4 = <IN4>;	
for (my $j = 0; $j < scalar @IN4; $j++){ # scan atom residue list (i.e. atom mask)
							my $IN3row = $IN3[$j];
							my @IN3row = split(/\s+/, $IN3row); 
			        my $IN4row = $IN4[$j];
							my @IN4row = split(/\s+/, $IN4row); 
							my $atom_numberQ = $IN3row[0];
							my $atom_labelQ = $IN3row[1];
							my $res_numberQ = $IN3row[2];
							my $res_labelQ = $IN3row[3];
							my $atom_numberR = $IN4row[0];
							my $atom_labelR = $IN4row[1];
							my $res_numberR = $IN4row[2];
							my $res_labelR = $IN4row[3];							
							  for (my $ll = 0; $ll < scalar @IN1; $ll++){ # scan corr files
								my $IN1row = $IN1[$ll];
								my @IN1row = split(/\s+/, $IN1row); 
			          my $res1 = $IN1row[0];
								my $atom1 = $IN1row[1];
								my $corr1 = $IN1row[2];
								my $IN2row = $IN2[$ll];
								my @IN2row = split(/\s+/, $IN2row); 
			          my $res2 = $IN2row[0];
								my $atom2 = $IN2row[1];
								my $corr2 = $IN2row[2];
								#if ($atom_numberQ == $atom1 && $corr1 =~/\d/){print "Q "."$res1\t"."$atom1\t"."$corr1\n";}
								#if ($atom_numberR == $atom2 && $corr2 =~/\d/){print "R "."$res2\t"."$atom2\t"."$corr2\n";}
								if ($atom_numberQ == $atom1 && $corr1 =~/\d/){print OUT1 "$res1\t"."$atom1\t"."$corr1\n";}
								if ($atom_numberR == $atom2 && $corr2 =~/\d/){print OUT2 "$res2\t"."$atom2\t"."$corr2\n";}
								}
							
							
							}					
	
	close IN1;
	close IN2;
	close IN3;
	close IN4;
}	



##############################################################################################################
print "\n\n collecting atomic correlation values (may take a minute)\n\n";
sleep(2);
open (OUT1, ">"."DROIDScorrelation.txt") or die "could not create output file\n";
print OUT1 "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."corr_ref\t"."corr_query\n";
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
							open(IN5, "<"."corr_mask_$fileIDr"."_$ii.txt") or die "could not open corr_mask_ file for $fileIDq\n";
              open(IN6, "<"."corr_mask_$fileIDq"."_$ii.txt") or die "could not open corr_mask_ file for $fileIDr\n";
	            my @IN5 = <IN5>;
              my @IN6 = <IN6>;
			        for (my $iii = 0; $iii < scalar @IN5; $iii++){
							    my $IN5row = $IN5[$iii];
									$IN5row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN5row = split(/\s+/, $IN5row);
									my $Qtest_atom_decimal = $IN5row[1];
									my $Qtest_atom = int($Qtest_atom_decimal);
							    #print "Q "."$Qtest_atom\t"."$atomnumberQ\n";
									if($atomnumberQ == $Qtest_atom){$corr_query = $IN5row[2];}
							  }	
							for (my $iii = 0; $iii < scalar @IN6; $iii++){
									my $IN6row = $IN6[$iii];
									$IN6row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN6row = split(/\s+/, $IN6row);
			            my $Rtest_atom_decimal = $IN6row[1];
									my $Rtest_atom = int($Rtest_atom_decimal);
									#print "R "."$Rtest_atom\t"."$atomnumberR\n";
									if($atomnumberR == $Rtest_atom){$corr_ref = $IN6row[2];}
							  }
							
					    if($resnumberR =~/\d/ && $corr_query=~/\d/ && $reslabelQ ne "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$corr_ref\t"."$corr_query\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$corr_ref\t"."$corr_query\n";
							    }
							if($resnumberR =~/\d/ && $reslabelQ eq "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$corr_ref\t"."NA\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$corr_ref\t"."NA\n";
							    }
							
							
							} 	
							
							close IN5;
              close IN6;
							
					
							}
					 
	
close IN3;
close IN4;
close OUT1;


###########################################################################################
print "\n\n parsing DROIDScorrelations by residue\n\n";
mkdir ("atomcorr") or die "please delete atomcorr folder from previous run\n";
open (IN, "<"."DROIDScorrelation.txt") or die "could not create output file\n";
my @IN = <IN>;
open (OUT2, ">"."DROIDScorrelationAVG.txt") or die "could not create output file\n";
print OUT2 "pos_ref\t"."res_ref\t"."res_query\t"."corr_ref_avg\t"."corr_query_avg\t"."delta_corr\t"."abs_delta_corr\n";
my @REFcorrAvg = ();
my @QUERYcorrAvg = ();
for (my $j = 0; $j < scalar @IN; $j++){ # scan atom type
			     my $INrow = $IN[$j];
	         my @INrow = split(/\s+/, $INrow); 
			     my $sample = $INrow[0];
					 my $pos_ref = $INrow[1];
					 my $res_ref = $INrow[2];
					 my $res_query = $INrow[3];
					 my $atomnumber = $INrow[4];
					 my $atomlabel = $INrow[5];
					 my $corr_ref = $INrow[6];
					 my $corr_query = $INrow[7];
					 push(@REFcorrAvg, $corr_ref);
					 push(@QUERYcorrAvg, $corr_query);
					 my $INnextrow = $IN[$j+1];
	         my @INnextrow = split(/\s+/, $INnextrow); 
			     my $next_pos = $INnextrow[1];
					 print OUT "$sample\t"."$pos_ref\t"."$res_ref\t"."$res_query\t"."$atomnumber\t"."$atomlabel\t"."$corr_ref\t"."$corr_query\n";
					 if ($homology eq "loose"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_query ne "xxx"){  # loose homology = collect all aligned residues  
           open (OUT, ">"."./atomcorr/DROIDScorrelation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."corr_ref\t"."corr_query\n";
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFcorrAvg);
					 $corr_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYcorrAvg);
					 $corr_query_avg = $statSCORE->mean();
					 $delta_corr = -($corr_ref_avg - $corr_query_avg);
					 $abs_delta_corr = abs($corr_ref_avg - $corr_query_avg);
					 if ($pos_ref =~ m/\d/ && $res_query ne "xxx"){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$corr_ref_avg\t"."$corr_query_avg\t"."$delta_corr\t"."$abs_delta_corr\n";}
					 if ($pos_ref =~ m/\d/ && $res_query eq "xxx"){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$corr_ref_avg\t"."NA\t"."NA\t"."NA\n";}
					 my @REFcorrAvg = ();
           my @QUERYcorrAvg = ();
					 if ($next_pos eq ''){next;}
					 }}
					 if ($homology eq "strict"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_ref eq $res_query && $res_query ne "xxx"){ # strict homology = collect only exact matching residues  
           open (OUT, ">"."./atomcorr/DROIDScorrelation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."corr_ref\t"."corr_query\n";
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@REFcorrAvg);
					 $corr_ref_avg = $statSCORE->mean();
					 $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
           $statSCORE->add_data (@QUERYcorrAvg);
					 $corr_query_avg = $statSCORE->mean();
					 $delta_corr = -($corr_ref_avg - $corr_query_avg);
					 $abs_delta_corr = abs($corr_ref_avg - $corr_query_avg);
					 if ($pos_ref =~ m/\d/ && $res_query ne "xxx"){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$corr_ref_avg\t"."$corr_query_avg\t"."$delta_corr\t"."$abs_delta_corr\n";}
					 if ($pos_ref =~ m/\d/ && $res_query eq "xxx"){print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$corr_ref_avg\t"."NA\t"."NA\t"."NA\n";}
					 my @REFcorrAvg = ();
           my @QUERYcorrAvg = ();
					 if ($next_pos eq ''){next;}
					 }}
																	
}
SKIP:
close IN;
close OUT;
sleep(2);
print "\n\n done parsing CPPTRAJ data files\n\n";
sleep(2);

system "perl GUI3_DROIDS.pl\n";	
}

##################################################################################################
