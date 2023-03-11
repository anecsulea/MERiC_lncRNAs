#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts unambiguously mapped reads.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################

## parameters 

my %parameters;
$parameters{"pathAllReads"}="NA";
$parameters{"pathUniqueReads"}="NA";


my @defaultpars=("pathAllReads","pathUniqueReads");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
	
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}


## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";


##############################################################

my $pathhits=$parameters{"pathAllReads"};
my $pathout=$parameters{"pathUniqueReads"};

##############################################################

## see if the file is zipped or not

my @ext=split("\\.",$pathhits);
my $nb=@ext;
my $extension=$ext[$nb-1];

### we take only the reads for which exactly one alignment was reported !!! 

my $input;

if($extension eq "gz"){
    open($input,"zcat $pathhits |");
}
else{
    if($extension eq "bam"){
	open($input,"samtools view -h $pathhits |");
    }
    else{
	open($input,$pathhits);
    }
}

open(my $output,">".$pathout);
my $line=<$input>;

my $nbkept=0;
my $nbtot=0;

while($line){
    
    chomp $line;

    my $firstchar=substr $line,0,1;

    if($firstchar eq "@"){
	print $output $line."\n";
	$line=<$input>;
	next;
    }
    
    my @s=split("\t",$line);
    $nbtot++;
    
    my $id=$s[0];
    my $nbhits="NA";

    for(my $u=11; $u<@s; $u++){
    	my @m=split(":",$s[$u]);
	if($m[0] eq "NH"){
	    $nbhits=$m[2]+0;
	    last;
	}
    }
    
    if($nbhits ne "NA"){
	if($nbhits == 1){
	    print $output $line."\n";
	    $nbkept++;
	}
    }
    else{
	print "could not find number of hits for ".$line."!!\n";
	exit(1);
    }

    $line=<$input>;
    
}

print "Kept ".$nbkept." reads with exactly one reported mapping out of ".$nbtot." total reads.\n";

close($input);
close($output);

##############################################################
