#!/usr/bin/perl -w

use strict;
use Spreadsheet::WriteExcel;
use File::Basename;
use Getopt::Std;

my @attributes = ("File","X2","X2red","PX2tot","XCG","CCG","SCG","XPh","CPh","SPh","XCh","CCh","SCh","XC","CC","SC","Xc1","Cc1","Sc1","Xc3","Cc3","Sc3","r","r12","RCG","RPh","Rm","sigR","A","nC2","nC1","VL","VHL","VCG","VPh","VCh","Vc2","Vc1","Vc3","D","DPP","DB","DC","DH1","dXH","dXH2","negP","dXH3");

my $dirName = $ARGV[0];
opendir(DIR, $dirName);
my @files = grep(/\.out$/,readdir(DIR));
closedir(DIR);

#Read an optional file name
shift @ARGV;
my %opt;
my $fileName;
#Look for a -f switch with an argument
getopt('f:', \%opt);
if(defined $opt{f}){
        $fileName = $opt{f};
        $fileName =~ /(^[^\.\/]+)/; #Ignore '.', '/', and anything that follows
	if(defined $1){
        	$fileName = $1 . ".xls";
        }
        #If no usable characters are extracted from the supplied file name,
        #use the default name fitdata.xls.
        else{
		$fileName = "fitdata.xls" 
	}
}
#If no file is given, name the output fitdata.xls
else{
	$fileName = 'fitdata.xls';
}
print "Using file name $fileName\n"; 
# Create a new Excel workbook
my $workbook = Spreadsheet::WriteExcel->new($fileName);
# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Add and define a format
my $format = $workbook->add_format(); # Add a format
$format->set_bold();
$format->set_align('center');


## Write attribute headers in the zeroth column
my $itemscol = 0;
my $itemsrow = 0;
foreach my $i (@attributes) {
    $worksheet->write($itemsrow++, $itemscol, $i, $format);
}
$itemscol++;

#The elements given by these indices of x(i) array will not be included in the output spreadsheet
#For details of what they correspond to, see an .out file: specifically, lines starting with set x(i)
#If a user decides to include more stuff in the spreadsheet, this list must be modified 
#according to his/her need.
my @ignoreThese = (27, 30, 31, 43, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 66);

foreach my $file (@files)
{
    my $itemsrow = 0;
    $file = $dirName."/".$file;
    open(OUT,"<$file");
    my $samplename = basename($file, ".out");
    print "wrote $samplename\n";
    $worksheet->write($itemsrow, $itemscol, $samplename);
    while(<OUT>) {    #This will read the file line by line
	    my $things = $_;

        #This block extracts X2, X2red, and PX2tot 
        if ($things =~ m/# MSD=/) {
            $things =~ m/X2=(.*) X2red/;
			$worksheet->write_number(++$itemsrow, $itemscol, $1);
			  
			$things =~ m/X2red=(.*) PT2/;
			$worksheet->write_number(++$itemsrow, $itemscol, $1);
			  
			$things =~ m/PX2tot=(.*) N/;
			$worksheet->write_number(++$itemsrow, $itemscol, $1);  
		}
        
        my $skip = 1;
        for (my $count = 0; $count < 68; $count++) {
            if ($things =~ m/set x\($count\) (.*);/) {
                #Check whether or not the line should be included in the output
                foreach my $i (@ignoreThese) {
                    if ($count == $i) {
                        $skip = 0;
                    }
                }     
                if ($skip) {
                    $worksheet->write_number(++$itemsrow,$itemscol,$1);
                }
            }
		}
	}
    $itemscol++;
    close(OUT);
}

#system("cp fitdata.xls $dirName");
#system("ls -l");
