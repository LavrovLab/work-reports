#!/usr/bin/perl
#  Extract the ends of sequences from fasta file.
#
# Version: 20060413
#
# Copyright 2006, Naoki Takebayashi <ffnt@uaf.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA

my $usage="\nUsage: $0 [-l length] [fasta_file [fasta_file ...]]\n\n".

    " -l: length of sequences ends you want. {default 30 bases]\n\n".

    "Reads in FASTA file(s) and extract the ends (connected by --).\n" .
    "If no lengths are specified, 30 bases are extracted by default,\n" .
    "If you want to extract 10 based from the beginning, and 20 from the end,\n".
    "specify \"-l 10,20\" (no space between the two integers).\n" .
    "If the sequence name contains tab, there will be a problem. Modify\n".
    "the value of \$sep in the script.\n";

my $sep = "\t";  # if you use tab in the sequence name, change this to
                 # other characters such as ","
my $replaceChar = "-"; # for -g, this character is used to replace

use Getopt::Std;
getopts('hl:') || die "$usage\n";

if (defined($opt_h)) {
    die "$usage\n";
}

die "$usage\n" if (@ARGV > 1);


my $endLen =30;
my $beginningLen=30;

if (defined($opt_l)) {
    if ($opt_l =~ /,/) {
	my @len = split (/,/, $opt_l);
	if (@len !=2) {
	    die "$usage\n";
	}
	if ($len[0] !~ /\d+/ || $len[1] !~ /\d+/) {
	    die "$usage\n";
	}
	($beginningLen, $endLen) = @len;
    } else {
	$beginningLen = $endLen = $opt_l;
    }
}


@ARGV = ('-') unless @ARGV; # take STDIN when no arg.

my $dnaFile;
while ($dnaFile = shift @ARGV)  {
# initialize the @seqArray, @seqNameArray, and $maxSeqLen
    my @dat = ReadInFASTA($dnaFile);
    my $numSeq = @dat;

    foreach $entry (@dat) {
	my @seq = split(/$sep/,$entry);
	print ">$seq[0]\n";
	my $head = substr($seq[1], 0, $beginningLen);
	my $tail = substr($seq[1], - $endLen);
	print "$head--$tail\n";
    }
}
exit(0);

#### functions

# takes an arg; name of a file from which data are read Then read in
# the data and make an array.  Each element of this array corresponds
# to a sequence: name tab data.
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
	    s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            s/[uU]/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . uc($_);
        }

	# checking no occurence of internal separator $sep.
	die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
	     "the input FASTA file contains this charcter. Make sure this " . 
	     "separator character is not used in your data file or modify " .
	     "variable \$sep in this script to some other character.\n")
	    if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
	$result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}
