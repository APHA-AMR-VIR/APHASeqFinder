#!/usr/bin/perl
#Badly constructed by Nick
# the file names MUST match perfectly or you get error!!!

## Must be run with python 3 environment
# Run from location of abricate_single output and seqfinder output
#	> perl ~/mnt/BactiPipes_MA2017/Scripts/abricate_seqfinder_combine/abricate_seqfinder_combine.pl 



system("echo 'Started running script at:'; date;");
system("echo 'Number of abricate files:'; ls *abricate | wc -l");
system("echo 'Number of SeqFinder files:'; ls *good_snps.csv | wc -l");
system("python /home/javi/APHASeqFinder/making_a_list_for_abricate_seqfinder_combine.py");
my $input = "./samples.csv"; 

open (LIST, "<$input") or die$!;

#place each line of the text file into an array.
my @list=<LIST>;
close LIST or die$!;
foreach my $line (@list) {
	chomp $line;
#print $line;
	my @words = split (",",$line);
	my @seq = split ("_",$line);
#create variables for file names.
	my $abricate = $words[0];
	my $seqfinder = $words[1];
#execute commands on the command line using "system command".
	print "Combining $abricate and $seqfinder \n";
	system("python /home/javi/APHASeqFinder/abricate_combine_with_seqfinder_v1.py $abricate $seqfinder");
}
system("echo 'Number of abricate files:'; ls *abricate | wc -l");
system("echo 'Number of SeqFinder files:'; ls *good_snps.csv | wc -l");
system("echo 'Number of joint files made:'; ls *abricate_seqfinder.csv | wc -l");	
system("python /home/javi/APHASeqFinder/compilation_for_abricate_seqfinder.py");	

end

