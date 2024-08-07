#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }

    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;

	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";

	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;

	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

#########################################################################################

sub reverseComplement{
    my $sequence=$_[0];

    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}

#########################################################################################

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
    }
}

###############################################################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $type=$s[2];

	if($type eq "exon"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $ss=$s[6];
	    my $strand="NA";

	    if($ss eq "+"){
		$strand="1";
	    } else{
		if($ss eq "-"){
		    $strand="-1";
		} else{
		    print "Weird strand!\n";
		    print $line."\n";
		    exit(1);
		}
	    }

	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $txid=findInfo("transcript_id", \@infoarray);

	    if($txid eq "NA"){
		print "weird! cannot find transcript_id in line.\n";
		print $line."\n";
		exit(1);
	    }

	    my $geneid=findInfo("gene_id", \@infoarray);

	    if($geneid eq "NA"){
		print "weird! cannot find gene_id in line.\n";
		print $line."\n";
		exit(1);
	    }

	    if(exists $transcripts->{$txid}){
		push(@{$transcripts->{$txid}{"start"}}, $start);
		push(@{$transcripts->{$txid}{"end"}}, $end);
	    } else{
		$transcripts->{$txid}={"start"=>[$start], "end"=>[$end], "chr"=>$chr, "strand"=>$strand, "gene"=>$geneid};
	    }
	}

	$line=<$input>;
    }

    close($input);
}

###############################################################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];

    my $res="NA";

    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    } else{
	my $nbreal=0;

	foreach my $g (@grepres){
	    $g =~ s/^\s+//; ## remove whitespace
	    my @u=split(" ",$g);

	    if($u[0] eq $pattern){
		$nbreal++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }

    if($nbreal>1){
	return "NA";
    }
    return $res;
}


##############################################################################

sub readGeneInfo{
    my $pathInfo=$_[0];
    my $info=$_[1];

    open(my $input,$pathInfo);

    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);

    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);

	my $id=$s[$header{"Gene stable ID"}];
	my $thisbio=$s[$header{"Gene type"}];

	$info->{$id}={"biotype"=>$thisbio};

	$line=<$input>;
    }

    close($input);
}

#############################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script extracts cDNA sequences in fasta format. \n";
    print "\n";
    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################################
##########################################################################################

my %parameters;

$parameters{"pathAnnotGTF"}="NA";
$parameters{"forbiddenChromo"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"forbiddenBiotypes"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathAnnotGTF", "forbiddenChromo", "pathGeneInfo", "forbiddenBiotypes", "pathGenomeSequence", "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];

    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;

	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
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

#####################################################################################
#####################################################################################

print "Reading transcripts...\n";

my %transcripts;
readGTF($parameters{"pathAnnotGTF"}, \%transcripts);
my $nbtx=keys %transcripts;
print "Found ".$nbtx." transcripts.\n";

print "Done.\n";

##############################################################

print "Reading gene info...\n";

my %geneinfo;
readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);
my $nbg=keys %geneinfo;

print "Found info for ".$nbg." genes.\n";

print "Done.\n";

##############################################################

my @fb=split(",", $parameters{"forbiddenBiotypes"});
my $nbf=@fb;

print "There are ".$nbf." forbidden biotypes: ".join(", ", @fb)."\n";

my %fbio;
foreach my $f (@fb){
    $fbio{$f}=1;
}

my %forbiddengenes;
my $nbg=0;

foreach my $gene (keys %geneinfo){
    my $bio=$geneinfo{$gene}{"biotype"};

    if(exists $fbio{$bio}){
	$forbiddengenes{$gene}=1;
	$nbg++;
    }
}

print "There are ".$nbg." genes with forbidden biotypes.\n";

##############################################################

print "Reading genome sequence...\n";
my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);
print "Done.\n";

##############################################################

print "Writing output...\n";

my @fc=split(",", $parameters{"forbiddenChromo"});
my $nbfc=@fc;

print "There are ".$nbfc." forbidden chromosomes: ".join(",", @fc)."\n";

my %forbidden;
foreach my $f (@fc){
    $forbidden{$f}=1;
}

open(my $output, ">".$parameters{"pathOutput"});

my $nbdiscarded=0;

foreach my $txid (keys %transcripts){
    my $chr=$transcripts{$txid}{"chr"};
    my $strand=$transcripts{$txid}{"strand"};
    my $geneid=$transcripts{$txid}{"gene"};

    if(exists $forbiddengenes{$geneid}){
	print "Excluding gene ".$geneid." ".$geneinfo{$geneid}{"biotype"}."\n";
	$nbdiscarded++;
    } else{
        if(!exists $forbidden{$chr}){
	    if(exists $genome{$chr}){
		my $nbexons=@{$transcripts{$txid}{"start"}};

		my %hashexons;

		for(my $i=0; $i<$nbexons; $i++){
		    my $start=${$transcripts{$txid}{"start"}}[$i];
		    my $end=${$transcripts{$txid}{"end"}}[$i];

		    if(exists $hashexons{$start}){
			print "Weird! already saw ".$start." for ".$txid."\n";
			exit(1);
		    }

		    $hashexons{$start}=$end;
		}

		my @uniquestart=keys %hashexons;
		my @sortedstart=sort {$a<=>$b} @uniquestart;

		my $sequence="";

		foreach my $start (@sortedstart){
		    my $end=$hashexons{$start};
		    my $thisseq=substr $genome{$chr}, ($start-1), ($end-$start+1);
		    $sequence.=$thisseq;
		}

		if($strand eq "-1" || $strand eq "-"){
		    $sequence=reverseComplement($sequence);
		}

		writeSequence($sequence, $geneid.":".$txid, $output);

	    } else{
		print "weird! cannot find ".$chr." in genome.\n";
	    }
	} else{
	    print "We discard ".$geneid." ".$txid.", wrong chromosome\n";
	    $nbdiscarded++;
	}
    }
}

print "We discarded ".$nbdiscarded." transcripts in total.\n";

close($output);

print "Done.\n";

##############################################################
