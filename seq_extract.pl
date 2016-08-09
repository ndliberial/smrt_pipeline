#!/usr/bin/perl
#This script performs 9 major tasks comprising read collection, mapping, filtering and assembling.
#Different procedures have been declared to handle each of these tasks which are activated on user click
#At the end of the assembly, the user now has the option of further processing or results using web based GUIs
 

use warnings;
use strict;
use Bio::SeqIO;
use LWP::Simple;
my $task = ""; #variable to hold task id to perform
my $ref_chr = ""; #reference chromosome declaration
my $cord1 = ""; #cordinate 1
my $cord2 = "";#cordinate 2
my $cnt = 0;#counter
my $custom_ref = ""; #name of file holding extracted sequence
print "\n!!!Welcome to the Analysis and Assembly Pipeline for the Recovery of Repetitive DNA using Single Molecule Sequence Reads!!!\n\n";
print "Please select from the options below your preferred task and then hit Enter.\n";
print "1. Start subread download\n2. Custom reference generation\n3. Create masked reference\n4. Align reads to masked reference\n5. Extract flanking reads from alignment\n6. Run pairwise alignment of flanking reads with custom reference\n7. Filter reads by score distribution\n8. Assemble reads\n9. Exit\n";

#accept task id
$task =  <STDIN>;
if ($task !~ m/[123456789]/i){
	print "Error! You can only select one of the listed tasks using the corresponding numbers!\n";
	exit;
}
else{
	if ($task == 1){
		&get_raw_reads;
	}
	elsif ($task == 2){
		&custom_extract;
	}
	elsif ($task == 3){
		&mask_ref;
	}
	elsif ($task == 4){
		&alignment;
	}
	elsif ($task == 5){
		&flank_reads;
	}
	elsif ($task == 6){
		&pairwise_alignment;
	}
	elsif ($task == 7){
		&filter_reads;
	}
	elsif ($task == 8){
		&assembly;
	}
	elsif ($task == 9){
		&finish;
	}
}
exit;
#subroutine to terminate pipeline
sub finish{
	print "Automated read processing and assembly now completed.\nTo determine the best reconstruction, ";
	exit;
}
#subroutine for performing assembly
sub assembly{
#this routine handles the assembly of reads initially filtered.
#user selects the distribution profile of choice and corresponding reads are then assembled with celera 8.2 assembler.
	print "Enter name/path of file containing reads to be assembled\n";
	my $choice = <STDIN>;
	chomp ($choice);
	unless (-e "$choice")   {print "Oops! Looks like file path does not exist.\nDo you wish to start again? Yes or No\n"; 
			$cnt =  <STDIN>;
			chomp($cnt);
			if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
				&assembly;
				return;
			}
			#otherwise terminate program
			else{
				exit;
			}
		}
		my @strings = split(/\./, $choice);
		$cnt = $strings[0];
		#generating formatted fragment file
		system "fastqToCA -libraryname chm1 -technology pacbio-raw -reads $choice > $cnt.frg";
		print "formatted fragment generation completed\nStarting assembly\n";
		system "runCA -d $cnt -p $cnt -s pac.spec $cnt.frg";
		print "Assembly finished. Check directory $cnt for details\n";
		print "Would you like to terminate this pipeline? Yes/No\n";
		$cnt = <STDIN>;
		chomp($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
		
			&assembly;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
		
	
}
#subroutine to handle read filtering
sub filter_reads{
	#Here, following pairwise alignment, a score distribution profile is applied to filter reads to include only reads with sufficient identity to ROI
	#This leads to the generation of 3 different read groups from the initial collection of flanking reads
	my $highest_scr = $_[1];
	my $output = "./last/last_scores.csv"; #alignment score file earlier generated
	my $flanked = $_[0]; #flanked reads
	my $aln_file = "./last/last_dist_mod.csv"; #read distribution profile
	my @score; #array to hold scores
	my $cnt = 0; #counter to keep track of read count
	my $rounded = 0; #variable to hold score rounded to significant values
	my $read98 = "./last/last98grp.fastq"; #98% read group containing fastq sequences
	my $read90 = "./last/last90grp.fastq"; #90% read group containing fastq sequences
	my $read100 = "./last/last100grp.fastq"; #100% read group containing fastq sequences
	my $read98txt = "./last/last98grp.txt"; #98% read group containing read ids
	my $read90txt = "./last/last90grp.txt"; #98% read group containing read ids
	my $read100txt = "./last/last100grp.txt"; #98% read group containing read ids
	#opening file handle to enable writing 
	open my $f98, '>>', $read98 or die "Could not open '$read98' $!\n";
	open my $f90, '>>', $read90 or die "Could not open '$read90' $!\n";
	open my $f100, '>>', $read100 or die "Could not open '$read100' $!\n";
	open my $f98txt, '>>', $read98txt or die "Could not open '$read98txt' $!\n";
	open my $f90txt, '>>', $read90txt or die "Could not open '$read90txt' $!\n";
	open my $f100txt, '>>', $read100txt or die "Could not open '$read100txt' $!\n";
	open my $n, '>>', $aln_file or die "Could not open '$aln_file' $!\n";
	my $scr = 0; #alignment score
	my $l98 = 0; #count of reads in the 98% group
	my $id = ""; #read id
	my $l90 = 0; #count of reads in the 90% group
	my $l100 = 0; #count of reads in the 100% group
	my $frequency = 0;
	if (!defined $flanked){
		print "Enter name/path of flanked reads\n";
		$flanked = <STDIN>;
		chomp ($flanked);
		print "Enter the highest score obtained from pairwise alignment\n";
		$highest_scr = <STDIN>;
		chomp ($highest_scr);
		unless (-e "$flanked")   {print "Oops! Looks like file path does not exist.\nDo you wish to start again? Yes or No\n"; 
			$cnt =  <STDIN>;
			chomp($cnt);
			if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
				&filter_reads;
				return;
			}
			#otherwise terminate program
			else{
				exit;
			}
		}
		#$flanked = $cnt;
		unless ($highest_scr =~ /^\d+$/) {
			print "Coordinates must only contain digits.\nDo you wish to start again? Yes or No\n";
			$cnt =  <STDIN>;
			chomp($cnt);
			
			if (($cnt eq "Yes") || ($cnt eq "y")){
				
				&filter_reads;
				return;
			}
			else{
				exit;
			}
		}
		#$highest_scr = $cnt;
		
	}
	$cnt = 0;
	print $n "\% Criteria\tScore distribution\tNo. of Reads\n"; #write header details of extraction
	for (my $i = 2; $i <= 100; $i+=2){ #set filtering criteria in even ratios
		
		$rounded = ($i/100) * $highest_scr; #determine score range corresponding to percentage of highest scoring alignment
		$rounded = int ( $rounded + 0.5 ); #round up percentage to nearest whole number
		print "$rounded\n"; #display score range for that percentile
		open my $n_, '<', $output or die "Could not open '$output' $!\n"; #open text file containing alignment scores
		while (<$n_>){
			$cnt++; #increment alignment score count
			@score = split (/\t/, $_); #split entry by tab
			$score[0] =~ s/>/@/; #replace fasta header with fastq header
			if ($score[1] =~ m/Score/){ #is the entry matching a heading? if so, do nothing
								
			}
			else{ #if not, then check that score is not greater than range 
				if (($score[1]> $scr) && ($score[1]<= $rounded)){
					$frequency++;
					
					
					if ($i > 2) #check that this is the 98 percentile group
					{
						open my $flank, '<', $flanked or die "Could not open '$flanked' $!\n"; #open flanking reads file
						while (my $j = <$flank>) { #start to read file read by read
							chomp $j;
							$score[0] =~ s/>/@/; #replace fasta header with fastq header
							if ($score[0] eq $j){#if read id matches entry, start filtering
								
								
								$j = <$flank>;
								$id = $score[0] . "\n" . $j;
								$j = <$flank>;
								$id = $id . $j;
								$j = <$flank>;
								$id = $id  . $j;
								#write fastq lines to corresponding file handle
								print $f98 "$id";
								$l98++; #increment count for 98% group
								last; #exit loop
								
								
							}
						}
						close ($flank); #close handle
						print $f98txt "$score[0]\t$score[1]"; #write read ids to corresponding file
					}
					if ($i > 10) #check that this is the 90 percentile 
					{
						
						open my $flank, '<', $flanked or die "Could not open '$flanked' $!\n"; #open flanked reads file for read
						while (my $j = <$flank>) { #loop through file for matching reads
							chomp $j; #remove space and new line characters
							
							$score[0] =~ s/>/@/; #replace fasta header with fastq header
							if ($score[0] eq $j){ #if read id matches read in flanked file, then begine extraction
								
								
								$j = <$flank>;
								$id = $score[0] . "\n" . $j;
								$j = <$flank>;
								$id = $id . $j;
								$j = <$flank>;
								$id = $id . $j;
								print $f90 "$id"; #write fastq entries to file handle
								$l90++; #increment count of read id
								last; #exit loop
								
							}
						}
						print $f90txt "$score[0]\t$score[1]"; # write read id to corresponding read group
						close ($flank); #close handle
					}
					if ($i > 0) #check that this is the 100 percentile
					{
						open my $flank, '<', $flanked or die "Could not open '$flanked' $!\n";  #open flanked reads file for read
						while (my $j = <$flank>) { #loop through file for matching reads
							chomp $j; #remove space and new line characters
							
							if ($score[0] eq $j){ #if read id matches read in flanked file, then begine extraction
								
								
								$j = <$flank>;
								$id = $score[0] . "\n" . $j;
								$j = <$flank>;
								$id = $id . $j;
								$j = <$flank>;
								$id = $id . $j;
								print $f100 "$id"; #write fastq entries to file handle
								$l100++; #increment count of read id
								last; #exit loop
								
							}
						}
						print $f100txt "$score[0]\t$score[1]"; # write read id to corresponding read group
						close ($flank); #close handle
					}
				
				}
			}
		}
		close ($n_); #close handle
		print $n "$i\t>$scr and <=$rounded\t$frequency\n"; #write details of distribution profile
		$scr = $rounded; #overwrite score 
		$frequency = 0; #initialise
	}
	print "90\% score distribution profile contains $l90 reads\n98\% score distribution profile contains $l98 reads\n100\% score distribution profile contains $l100 reads\n"; #display to screen a summary of filtering
	print "Would you like to proceed with assembly? Yes/No\n";
	$cnt = <STDIN>;
	chomp($cnt);
	if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
		
			&assembly;
			return;
	}
		#otherwise terminate program
	else{
		exit;
	}
}
#subroutine for pairwise alignment and score generation
sub pairwise_alignment{
#Using LAST 4.75 aligner, this routine performs pairwise alignment of each extracted read with the custom reference and stores the reads with alignment score > 0

	my $search_str_file = ""; #variable to hold temporary fasta sequence
	my $query_file = $_[0]; #extracted read file
	my $search_str = 0;
	my $output = "./last/last_scores.csv"; #output path for alignment scores
	my $aln_file = "./last/t_align.fasta"; #temporary file for reads
	my @score; #array to hold alignment scores
	my $id = ""; #variable to hold read id
	my $cnt = 0; #counter to keep track of read count
	my $final_entry = ""; #variable to hold aligned read id
	my $flag = 0; #variable to keep track of read selection
	my $lcnt = 0; #variable to track unaligned reads
	my $scr = 0; #alignment score
	my $c = 0; #variable to track aligned reads
	my $hscr = 0; #highest alignment score
	my $dbref = "./"; #indexed masked reference
	#$query_file = $_[0];
	if (!defined $query_file){
		print "Enter name/path of flanked reads\n";
		$query_file = <STDIN>;
		chomp ($query_file);
		unless (-e "$query_file")   {print "Oops! Looks like file path does not exist.\nDo you wish to start again? Yes or No\n"; 
			$cnt = <STDIN>;
			chomp($cnt);
			if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
				
				&pairwise_alignment;
				return;
			}
			#otherwise terminate program
			else{
				exit;
			}
		}
	}
	open my $fh1, '<', $query_file or die "Could not open '$query_file' $!\n"; #file handle for extracted read file
	open my $n_, '>>', $output or die "Could not open '$output' $!\n"; #file handle for writing alignment output
	
	print "Enter name/path of custom reference\n";
	$cnt = <STDIN>;
	chomp ($cnt);										
	unless (-e "$cnt")   {print "Oops! Looks like file path does not exist.\nDo you wish to start again? Yes or No\n"; 
		chomp ($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&pairwise_alignment;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	}
	system "lastdb $dbref $cnt"; #using last to index custom reference
	$cnt = 0;
	print $n_ "Read_ID\tAlignment Score\n";
	 while (my $l = <$fh1>) { #start reading extracted read file
		 chomp $l; #remove white space or new line character
		if ($l =~ /^@(\S+)/){ #check that line is a valid fastq header
			$cnt++; #increment read count
			$flag = 1; #set flag to 1 as fastq header found
			$id = $l; #assign read id
			$id =~ s/@/>/; #replace fastq header with fasta header
			$search_str_file =  $cnt . ".fasta"; #file name to hold each pairwise alignment result
		}
		elsif ($l =~ m/[ATCG]/){ #check for valid fasta sequence
			if ($flag == 1){ #confirm that read sequence matches previously identified read id
				
				open my $temp, '>', $aln_file or die "Could not open '$aln_file' $!\n"; #open temporary fasta file 			
				print $temp "$id\n$l"; #write fasta sequence to temporary file
				system "lastal -k1 -T0 -m10 -w0 -g1.0 $dbref $aln_file > $search_str_file"; #perform last alignment and store output
				$final_entry = $id; #reassign read id
				$final_entry =~ s/>/@/; #replace fasta header with fastq header
				open my $fh, '<', $search_str_file or die "Could not open '$search_str_file' $!\n"; #open alignment file 
				 while (my $line = <$fh>) { #start reading line by line
					chomp $line; #remove space or new line character
					if ($line =~ /score/){ #is this the line with alignment score entry?
						@score = split(/=/, $line); #if so, separate score from other entries
						if ($score[1] > $scr){ #is this score greater than previous score? 
							$scr = $score[1]; #if so, overwrite previous score
							
							
						}
						
					}
					#if line dose not match alignment score entry, ignore
					else{
						
					}
				}
				if ($scr > 0){ #if score is greater than zero, store read id and score value in file
					print $n_ "$final_entry\t$scr\n";
					print "highest score in file $search_str_file is $scr\n"; #print to screen the highest score per aligned read 
					if ($scr > $hscr){ #check if next read has a higher score than previous read
						$hscr = $scr; #if so, replace score
					}
					$scr = 0; #initialise alignment score
					$c++; #increment count of aligned reads
				}
				else{
					$lcnt++; #increment count of unaligned read
				}	
			
				
				close ($fh); #close file handle for alignment file
				close ($temp); #close file handle for temporary fasta file
			}
			#system "rm $search_str_file"; #remove temporary fasta file
		}
	 
		else{
		
		}
		
		
	 }
	 #display to screen a summary report on the pairwise alignment
	 print "Total number of initial reads is $cnt\nTotal number of unaligned reads is $lcnt\nTotal number of aligned reads is $c\n";
	 print "The highest final score is $hscr\n";
	
}
#subroutine for read mapping
sub alignment{
	#This routine gives the option to either perform read mapping locally or on the cluster. 
	#Mapping locally runs simply without any further work but on the cluster requires the use of a pbs compliant script.
	#Depending on your cluster machine and available resources, you will need to modify the attached sample job schedule script
	print "1. To run alignment locally, select 1\n2. To run alignment on a cluster, select 2\n";
	$cnt = <STDIN>;
	chomp ($cnt);
	if ($cnt == 1){
		print "Enter name/path of masked reference\n";
		$cnt = <STDIN>;
		chomp($cnt);
		unless (-e "$cnt")   {print "Oops! Looks like file path does not exist.\nDo you wish to start again? Yes or No\n"; 
			if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
				&alignment;
				return;
			}
			#otherwise terminate program
			else{
				exit;
			}
		}
		print "Indexing masked reference started...\n";
		system "bwa index $cnt";
		print "Indexing completed\n";
		#system "cd datasets/";
		#system "pwd";
		for (my $i = 1; $i <= 20; ++$i){
			print "Mapping subread $i with reference\n...";
			system "bwa mem -x pacbio $cnt ./raw_seq/$i.fq > $i.sam";
			print "Mapping completed\nNow converting, sorting and indexing output to bam\n";
			system "samtools view -bS $i.sam > $i.bam";
			system "samtools sort -m 4000M $i.bam $i.sorted";
			system "samtools index $i.sorted.bam";
			system "mv *.bam *.bai ./datasets";
			system "rm *.sam";
			print "Conversion completed\n";
			
			
		}
		
	}
	elsif ($cnt ==2){
		
	}
	else{
		print "Error! You can only choose one of the listed tasks using the corresponding numbers!\n";
		exit;
	}
	print "Would you like to proceed with flanking read extraction? Yes/No\n";
	$cnt = <STDIN>;
	chomp($cnt);
	if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
		
			&flank_reads;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
}
#subroutine to automatically download all pacbio (P5_C3) CHM1 htert reads 
sub get_raw_reads{
#This perl script automates the process of downloading raw SMRT reads sequenced to 54X coverage from PacBio's online data repository using wget
	print "Supply full path to text file containing download links\n";
	my $filelinks = <STDIN>; #get text file containing links to online data repository for download
	chomp($filelinks);#removing new line and space characters from file name
	unless (-e "$filelinks")   {print "Oops! Looks like the file path does not exist.\nDo you wish to start again? Yes or No\n"; 
		$cnt =  <STDIN>;
		chomp($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&get_raw_reads;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	}
	my $cnt = 0 ; #counter to track number of subreads being downloaded
	my $f_path = $filelinks; #store text file containing download links to another variable
	my $store_path = "./raw_seq/"; #specify directory for storing reads downloaded
	my $file_name = "";#variable to hold downloaded read 
	
	#open file containing download links for reading
	open my $fh, '<', $filelinks or die "Could not open '$filelinks' $!\n";
		 
	#loop through links file and download each subread
	
	while (my $line = <$fh>) {
		chomp $line; #trim off all trailing white spaces
		$cnt++; #increment count of files downloaded
		$file_name = $store_path . $cnt . ".fq"; #set full storage path for each subread using the value of count 
		print "Download started...\n";
		system "wget -q $line -O $file_name"; #using system call to initiate the download process via wget command
		print "Subread number $cnt successfully downloaded\n";
		
	}
	print "Subread download finished\n";
	print "Would you like to proceed with custom reference generation? Yes/No\n";
	$cnt = <STDIN>;
	chomp($cnt);
	if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&custom_extract;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	
}
#subroutine to extract 2.5kb flanking reads from alignment
#for this subroutine to work, both sorted bam and index files have to be present
sub flank_reads{
	my $sorted_bam = ""; #holds the sorted bam from which sequences are to be extracted
	my $raw_extract = ""; #palceholder for initial extraction
	my $final_extract = ""; #variable to hold filtered read id
	my $uniq_final_extract = ""; #variable to hold unique read ids
	print "Supply the start position of interest. This ideally should be the start coordinate of the ROI\n";
	my $start_pos = <STDIN>;
	chomp($start_pos);
	$start_pos -= 2500; #flanking position from the start
	print "Supply the end position of interest. This ideally should be the end coordinate of the ROI\n";
	my $end_pos = <STDIN>;
	chomp($end_pos);
	 $end_pos += 2500; #flanking position from the end
	my @strings; #array to hold splitted strings from initial bam extract by lines
	print "Give the chromosome ID of interest\n";
	my $chr = <STDIN>; #variable to hold chromosome
	chomp($chr); 
	#looping through all sorted alignments to extract reads mapped to region of interest
	for (my $i = 1; $i <= 18; $i+=1){
		$raw_extract = "./datasets/" . $i . "initial_flank.fa";
		$final_extract = "./datasets/" .  $i . "_mapped_flanking.fa";
		$uniq_final_extract =  "./datasets/" . $i . "uniq_mapped_flanking.fasta";
		$sorted_bam =  "./datasets/" . $i . '.sorted.bam';
			
		#extracting read id from sorted bam
		system "samtools view $sorted_bam $chr:$start_pos-$end_pos | cut -f 1 > $raw_extract";
		
		#reading through read extract and formatting 
		open (IN, "$raw_extract");
		open my $outfile, '>>', $final_extract or die "Could not open '$final_extract' $!\n";
		
		while (<IN>){
			chomp $_;
			@strings = split(/\t/, $_);
			print $outfile "\>$strings[0]\n";
		}
		
		#using sort function to remove duplicate read ids
		system "sort $final_extract | uniq > $uniq_final_extract";
		#removing temporary files from directory
		system "rm $raw_extract";
		system "rm $final_extract";
		#print "values are $raw_extract, $final_extract and $sorted_bam\n"; exit;
		&fastq_extract($uniq_final_extract, $i);
			 
	}
	system "cat ./flank/*.fq > ./flank/flank.fastq"; #merge all extracted files into a single fastq file
	system "rm ./flank/*.fq"; #delete individual fastq files
	system "rm ./datasets/*.fasta"; #delete read id list
	print "Would you like to proceed with pairwise alignment? Yes/No\n";
	$cnt = <STDIN>;
	chomp($cnt);
	if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
		
			&pairwise_alignment("flank.fastq");
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	return;
}
#subroutine for extracting flanking read sequence from raw fastq file
sub fastq_extract{
	#taking the flanking read ids from the sburoutine flank_reads, corresponding fastq formatted sequences are extracted from the original downloaded fastq files
	my $readidfilename = $_[0]; #specify read id list file as argument 
	my $fastqfilename; #downloaded fastq file
	my $file1 =  $_[1]; #output location
	my %readids; #perl hash to store read ids
	my $c = 0; #variable to count number of read ids processed
	my $cnt = 0; #variable to count number of reads extracted
	my @strings; #perl array to formatted read id
	
	$fastqfilename = "./raw_seq/$file1.fq"; #path to original fastq file
	$file1 = "./flank/$file1.fq";#location to store extracted reads
	open( my $readidfile, "<", $readidfilename );
	open( my $fastqfile, "<", $fastqfilename );
	
	open my $output, '>>', $file1 or die "Could not open '' $!\n";
	while (<$readidfile>) {
		
		chomp();
		$_  =~ s/>/@/; #if read id is preceeded by a fasta identifier, replace with fastq identifier
		
		
		$readids{$_} = ''; #store read id in a hash
		
	}
	#loop through original fastq sequence to find read id match
	while (<$fastqfile>) {
		
		chomp();
		my $line = $_;
		if ($line =~ /^@/ ) { #check that line matches standard fastq read id line
	
			@strings = split(/ /, $line); #if so, split line
			$line = $strings[0]; #store read id
		
			chomp($line);
			if (exists( $readids{$line} )) { #check if read id exist in hash already
				#extract all four lines of standard fastq read and write to extracted file
				print {$output}( $line . "\n" );
				$line = <$fastqfile>;
				print {$output}($line);
				$line = <$fastqfile>;
				print {$output}($line);
				$line = <$fastqfile>;
				print {$output}($line);	
				$cnt++;
				print "fasta file is $readidfilename and fastq file $fastqfilename and $file1\n";
				
				
			}
			#imcrement count of matching read id
			$c++;
		}
	
	}
	print "Original fastq file contained $c number of reads\n";
	print "Flanking extract contains $cnt number of reads\n";
	$cnt = 0;
	return;
	
		
}
#subroutine to mask reference using bedtools
sub mask_ref{
	print "Enter full path of reference to be masked\n";
	my $read_file =  <STDIN>; #unmasked reference 
	print "Enter output path for the masked reference\n";
	my $write_file = <STDIN>; #masked reference
	print "Enter full path of the bed formatted text file containing coordinates\n";
	my $masked_region =  <STDIN>; #bed file containing cordinates to be masked
	chomp($read_file); chomp($write_file); chomp($masked_region); #removing white space and newline characters from file names
	unless ((-e "$read_file")  && (-e "$masked_region"))  {print "Oops! Looks like a file path does not exist.\nDo you wish to start again? Yes or No\n"; 
	$cnt =  <STDIN>;
	chomp($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&mask_ref;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	}
	#using system call to initiate reference masking with the help of bedtools
	system "bedtools maskfasta -fi $read_file -fo $write_file -bed $masked_region";
	print "Reference successfully masked\n";
	print "Would you like to proceed with read mapping? Yes/No\n";
	$cnt = <STDIN>;
	chomp($cnt);
	if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&alignment;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	return;
}
#this subroutine returns dna sequence matching region specified by coordinates for the gievn reference chromosome
#four arguments are required for this subroutine to run perfectly and they are; start and stop coordinates, full path of reference chromosome, and path for output
sub custom_extract{
	my $input = ""; #sequence line	
	print "Enter the start coordinate\n";
	$cord1 = <STDIN>;
	print "Enter the stop coordinate\n";
	$cord2 = <STDIN>;
	chomp($cord1); chomp($cord2); #removing white space and newline characters 
	#check that coordniates match only numbers
	unless (($cord1 =~ /^\d+$/) && ($cord2 =~ /^\d+$/)) {
		print "Coordinates must only contain digits.\nDo you wish to start again? Yes or No\n";
		$cnt =  <STDIN>;
		chomp($cnt);
		
		if (($cnt eq "Yes") || ($cnt eq "y")){
			#print "count is $cnt\n";
			&custom_extract;
			return;
		}
		else{
			exit;
		}
	}
	print "Enter full path of reference chromosome\n";
	$ref_chr = <STDIN>;
	chomp($ref_chr); #removing white space and newline characters from file names
	
	
	unless (-e "$ref_chr") {print "Oops! File path does not exist.\nDo you wish to start again? Yes or No\n"; 
	$cnt =  <STDIN>;
	chomp($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
			&custom_extract;
			return;
		}
		#otherwise terminate program
		else{
			exit;
		}
	}
	# open reference for read
	open(IN,"$ref_chr") or die "can't open $ref_chr"; 
	while (<IN>) { 
		if ($_ =~ /^>(\S+)/){
		}
		else{
			$input .= $_;
		}
	}
	close (IN);
	print "Enter name/path for output file\n";
	$custom_ref = <STDIN>;
	chomp($custom_ref); #removing white space and newline characters from file names
	#Using Bio::SeqIO to perform sequence extraction
	my $in = Bio::SeqIO->new(-file => "$ref_chr",
                         -format => 'Fasta');
    my $seq = $in->next_seq();
   
   # print "end is $cord2, start is $cord1 and length is $input\n"; exit;
    if (($cord2 > $cord1) && (length($input) > $cord2)) {
    	#Extract a subsequence
    	my $subseq = $seq->subseq($cord1, $cord2);
    	#opening file for output
    	open my $f, '>>',  "$custom_ref" or die "Could not open '' $!\n";
    	#printing extracted sequence to file
    	print $f ">custom:$cord1-$cord2\n$subseq\n";
    	#printing length of extraction to screen
    	print "Length of extracted region is " . length($subseq) . "\n";
    	print "Would you like to proceed with reference masking? Yes/No\n";
		$cnt = <STDIN>;
		chomp($cnt);
		if (($cnt eq "Yes") || ($cnt eq "y")){ #if user chooses to repeat the task, then start afresh
				&mask_ref;
				return;
			}
			#otherwise terminate program
			else{
				exit;
			}
    	return;
    }
    else{
    	print "Sorry! Kindly check that end position is not less than start position or greater than sequence.\n";
    	return;
    }

}