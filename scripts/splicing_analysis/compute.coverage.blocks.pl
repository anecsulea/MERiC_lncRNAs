#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readExonBlocks{
    my $pathExons=$_[0];
    my $exons=$_[1];
    my $geneex=$_[2];

    open(my $input,$pathExons);
    
    my $line=<$input>; ## header
       
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $gene=$s[0];

	my $chr=$s[2];
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	my $id=$chr.",".$start.",".$end.",".$strand;
	
	$exons->{$id}={"chr"=>$chr,"start"=>$start,"end"=>$end,"strand"=>$strand,"gene"=>$gene};

	if(exists $geneex->{$gene}){
	    push(@{$geneex->{$gene}},$id);
	}
	else{
	    $geneex->{$gene}=[$id];
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub orderExons{
    my $exons=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $strand=$exons->{$exid}{"strand"};
	my $start=$exons->{$exid}{"start"};
	my $end=$exons->{$exid}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$exid);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		push(@{$hashpos{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};

	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################

sub computeCoverageExons{
    my $exons=$_[0]; ## ordered coordinates
    my $coverage=$_[1]; ## ordered coordinates
    my $strand=$_[2];
    my $covexons=$_[3];

    foreach my $chr (keys %{$exons}){
	my $nbex1=@{$exons->{$chr}{"start"}};
	
	if(exists $coverage->{$chr}){
	    my $nbex2=@{$coverage->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$exons->{$chr}{"start"}}[$i];
		my $end1=${$exons->{$chr}{"end"}}[$i];
		my $strand1=${$exons->{$chr}{"strand"}}[$i];
		my $id1=${$exons->{$chr}{"id"}}[$i];
		
		my $length=$end1-$start1+1.0;

		my $j=$firstj;
		
		while($j<$nbex2 && ${$coverage->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		my $sumcov=0;
		my $covlen=0;

		$firstj=$j;
		
		if(($strand1 eq $strand)){

		    while($j<$nbex2 && ${$coverage->{$chr}{"start"}}[$j]<=$end1){
			
			my $start2=${$coverage->{$chr}{"start"}}[$j];
			my $end2=${$coverage->{$chr}{"end"}}[$j];
			my $thiscov=${$coverage->{$chr}{"cov"}}[$j];
			
			my $M=max($start1,$start2);
			my $m=min($end1,$end2);
			
			if($M<=$m){
			    $sumcov+=($thiscov+0.0)*($m-$M+1);
			    $covlen+=($m-$M+1);
			}

			$j++;
		    }
		    
		    my $meancov=$sumcov/$length;

		    $covexons->{$id1}={"meancov"=>$meancov,"coveredlength"=>$covlen};
		}
	    }
	}
	else{
	    for(my $i=0; $i<$nbex1; $i++){
		my $id1=${$exons->{$chr}{"id"}}[$i];
		my $strand1=${$exons->{$chr}{"strand"}}[$i];
		
		if($strand1 eq $strand){
		    $covexons->{$id1}={"meancov"=>0,"coveredlength"=>0};
		}
	    }
	}
    }
}

##############################################################

sub readCoverage{

    my $path=$_[0];
    my $refcov=$_[1];
    
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    } else{
	open($input,$path);
    }

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);
	my $chr=$s[0];
	my $start=$s[1]+1;
	my $end=$s[2];
	my $cov=$s[3]+0;

	if($cov>0){
	    if(exists $refcov->{$chr}){
		push(@{$refcov->{$chr}{"start"}},$start);
		push(@{$refcov->{$chr}{"end"}},$end);
		push(@{$refcov->{$chr}{"cov"}},$cov);
	    }
	    else{
		print $chr."\n";
		$refcov->{$chr}={"start"=>[$start],"end"=>[$end],"cov"=>[$cov]};
	    }
	}

	$line=<$input>;
    }
    
    close($input);
    
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes coverage for exons. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathCoverageForward"}="NA";
$parameters{"pathCoverageReverse"}="NA";
$parameters{"pathOutputExons"}="NA";
$parameters{"pathOutputGenes"}="NA";


my %defaultvalues;
my @defaultpars=("pathExonBlocks","pathCoverageForward","pathCoverageReverse","pathOutputExons","pathOutputGenes");

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

#####################################################################
#####################################################################

print "Reading exon blocks...\n";
my %exons;
my %genes;
readExonBlocks($parameters{"pathExonBlocks"},\%exons,\%genes);
print "Done.\n";

print "Ordering exons...\n";

my %orderedexons;
orderExons(\%exons, \%orderedexons);

print "Done.\n";

##############################################################

my %coveragesense;
my %coverageantisense;

print "Forward strand\n";

print "Reading coverage...\n";

my %coverageforward;
readCoverage($parameters{"pathCoverageForward"},\%coverageforward);
print "Done.\n";

print "Computing coverage for exons...\n";
computeCoverageExons(\%orderedexons, \%coverageforward, "1",\%coveragesense);
computeCoverageExons(\%orderedexons, \%coverageforward, "-1",\%coverageantisense);
print "Done.\n";

%coverageforward=();

print "Reverse strand\n";

print "Reading coverage...\n";
my %coveragereverse;
readCoverage($parameters{"pathCoverageReverse"},\%coveragereverse);
print "Done.\n";

print "Computing coverage for exons...\n";
computeCoverageExons(\%orderedexons, \%coveragereverse, "-1",\%coveragesense);
computeCoverageExons(\%orderedexons, \%coveragereverse, "1",\%coverageantisense);
print "Done.\n";

%coveragereverse=();

print "Done.\n";

##############################################################

print "Witing output...\n";

my $pathex=$parameters{"pathOutputExons"};

my $output;

if($pathex ne "NA"){
    open($output,">".$parameters{"pathOutputExons"});
}
open(my $outputgene,">".$parameters{"pathOutputGenes"});

print $outputgene "GeneID\tMeanCoverageSense\tCoveredLengthSense\tMeanCoverageAntisense\tCoveredLengthAntisense\n";

my %covgenes;

foreach my $gene (keys %genes){
    my $exoniclength=0;

    my $sumcovsense=0;
    my $totcovlensense=0;

    my $sumcovantisense=0;
    my $totcovlenantisense=0;
    
    foreach my $idex (@{$genes{$gene}}){
	my $chr=$exons{$idex}{"chr"};
	my $strand=$exons{$idex}{"strand"};
	my $start=$exons{$idex}{"start"};
	my $end=$exons{$idex}{"end"};

	$exoniclength+=($end-$start+1);

	my $meancovsense=0;
	my $covlensense=0;

	if(exists $coveragesense{$idex}){
	    $meancovsense=$coveragesense{$idex}{"meancov"};
	    $covlensense=$coveragesense{$idex}{"coveredlength"};
	}

	$sumcovsense+=($meancovsense+0.0)*($end-$start+1.0);
	$totcovlensense+=$covlensense;
	

	my $meancovantisense=0;
	my $covlenantisense=0;

	if(exists $coverageantisense{$idex}){
	    $meancovantisense=$coverageantisense{$idex}{"meancov"};
	    $covlenantisense=$coverageantisense{$idex}{"coveredlength"};
	}

	$sumcovantisense+=($meancovantisense+0.0)*($end-$start+1.0);
	$totcovlenantisense+=$covlenantisense;

	if($pathex ne "NA"){
	    print $output $gene."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$meancovsense."\t".$covlensense."\t".$meancovantisense."\t".$covlenantisense."\n";
	}
	
    }

    my $meancovgenesense=($sumcovsense+0.0)/($exoniclength+0.0);
    my $meancovgeneantisense=($sumcovantisense+0.0)/($exoniclength+0.0);

    print $outputgene $gene."\t".$meancovgenesense."\t".$totcovlensense."\t".$meancovgeneantisense."\t".$totcovlenantisense."\n";
}

if($pathex ne "NA"){
    close($output);
}

close($outputgene);

print "Done.\n";

##############################################################
