#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

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

#############################################################################################

sub extractPromoterCoordinates{
    my $txcoords=$_[0];
    my $winsize=$_[1];
    my $promotercoords=$_[2];

    foreach my $txid (keys %{$txcoords}){
	my $chr=$txcoords->{$txid}{"chr"};
	my $strand=$txcoords->{$txid}{"strand"};
	my $txstart=min @{$txcoords->{$txid}{"start"}};
	my $txend=max @{$txcoords->{$txid}{"end"}};

	my $promstart;
	my $promend;
	
	if($strand eq "+" || $strand eq "1"){
	    $promstart=$txstart-$winsize;
	    $promend=$txstart+$winsize;
	} else{
	    if($strand eq "-" || $strand eq "-1"){
		$promstart=$txend-$winsize;
		$promend=$txend+$winsize;
	    } else{
		print "Weird! don't know what to do for strand ".$strand."\n";
		exit(1);
	    }
	}
    
	$promotercoords->{$txid}={"chr"=>$chr, "start"=>$promstart, "end"=>$promend};
    }
}

#############################################################################################

sub readPCHiCContacts{
    my $pathin=$_[0];
    my $unorderedfrag=$_[1];
    my $type=$_[2];
    my $contacts=$_[3];
    
    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	
	my @s=split("\t", $line);

	my $thistype=$s[$header{"type"}];

	if($type eq "all" || $thistype eq $type){
	    ## bait
	    my $chrbait=$s[$header{"chr_bait"}];
	    my $prefixbait=substr $chrbait, 0, 3;
	    
	    if($prefixbait eq "chr"){
		$chrbait=substr $chrbait, 3;
	    }
	    
	    my $startbait=$s[$header{"start_bait"}]+0;
	    my $endbait=$s[$header{"end_bait"}]+0;

	    my $idbait=$chrbait.":".$startbait.":".$endbait;

	    if(!exists $unorderedfrag->{$idbait}){
		$unorderedfrag->{$idbait}={"chr"=>$chrbait, "start"=>$startbait, "end"=>$endbait};
	    }

	    ## other fragment
	    my $chrother=$s[$header{"chr"}];
	    my $prefixother=substr $chrother, 0, 3;
	    
	    if($prefixother eq "chr"){
		$chrother=substr $chrother, 3;
	    }
	    my $startother=$s[$header{"start"}]+0;
	    my $endother=$s[$header{"end"}]+0;
	    
	    my $idother=$chrother.":".$startother.":".$endother;

	    if(!exists $unorderedfrag->{$idother}){
		$unorderedfrag->{$idother}={"chr"=>$chrother, "start"=>$startother, "end"=>$endother};
	    }

	    ## store contact
	    if(exists $contacts->{$idbait}){
		$contacts->{$idbait}{$idother}=1; 
	    } else{
		$contacts->{$idbait}={$idother=>1}; 
	    }
	}
	
	$line=<$input>;
    }

    close($input);
}

#############################################################################################

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

##################################################################################################

sub orderCoords{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $id (keys %{$coords}){
	my $chr=$coords->{$id}{"chr"};
	my $start=$coords->{$id}{"start"};
	my $end=$coords->{$id}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$id);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$id],"end"=>[$end]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$id],"end"=>[$end]}};
	}
    }
   

    foreach my $chr (keys %hashpos){
	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];
		
		if(!(exists $ordered->{$chr})){
		    $ordered->{$chr}={"start"=>[],"end"=>[],"id"=>[]};
		}

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##################################################################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $overlap=$_[2];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nb1=@{$coords1->{$chr}{"start"}};
	    my $nb2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nb1; $i++){
		
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];
			
		my $j=$firstj;
		
		while($j<$nb2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nb2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];
		    
		    my $M=max($start1,$start2);
		    my $m=min($end1,$end2);
		    
		    if($M<=$m){
			my $fr1=($m-$M+1)/($end1-$start1+1);
			my $fr2=($m-$M+1)/($end2-$start2+1);
			
			if(exists $overlap->{$id1}){
			    $overlap->{$id1}{$id2}=1;
			}
			else{
			    $overlap->{$id1}={$id2=>1};
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}

##################################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script extracts promoter-promoter contacts from PCHi-C data. \n";
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
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathInteractions"}="NA";
$parameters{"windowSize"}="NA";
$parameters{"interactionType"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathAnnotGTF", "pathGeneInfo", "pathInteractions", "windowSize", "interactionType", "pathOutput");

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

print "Reading transcript coordinates...\n";

my %txcoords;

readGTF($parameters{"pathAnnotGTF"}, \%txcoords);

print "Done.\n";

print "Reading gene info...\n";

my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

print "Done.\n";

#####################################################################################

print "Extracting promoter coordinates...\n";

my $winsize=$parameters{"windowSize"}+0;

print "window size: ".$winsize."\n";

my %promotercoords;

extractPromoterCoordinates(\%txcoords, $winsize, \%promotercoords);
 
print "Done.\n";

#####################################################################################

print "Reading PCHi-C contacts...\n";

my %allfragments;
my %contacts;

my $type=$parameters{"interactionType"};

print "Keeping ".$type." interactions.\n";

readPCHiCContacts($parameters{"pathInteractions"}, \%allfragments, $type, \%contacts);

my $nbfrag=keys %allfragments;

print "Saw ".$nbfrag." fragments in total.\n";

print "Done.\n";

#####################################################################################

print "Ordering promoter coordinates...\n";

my %orderedpromoters;

orderCoords(\%promotercoords, \%orderedpromoters);

print "Done.\n";

print "Ordering fragment coordinates...\n";

my %orderedfragments;

orderCoords(\%allfragments, \%orderedfragments);

print "Done.\n";

#####################################################################################

print "Extracting overlap between restriction fragments and promoters...\n";

my %overlap;

extractOverlap(\%orderedfragments, \%orderedpromoters, \%overlap);

print "Done.\n";

#####################################################################################

print "Analyzing overlap and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID1\tGeneID1\tGeneType1\tChr1\tPromoterStart1\tPromoterEnd1\tStrand1\tTranscriptID2\tGeneID2\tGeneType2\tChr2\tPromoterStart2\tPromoterEnd2\tStrand2\n";

foreach my $idbait (keys %contacts){
    if(exists $overlap{$idbait}){
	foreach my $idother (keys %{$contacts{$idbait}}){
	    if(exists $overlap{$idother}){

		foreach my $idtx1 (keys %{$overlap{$idbait}}){
		    my $gene1=$txcoords{$idtx1}{"gene"};
		    my $strand1=$txcoords{$idtx1}{"strand"};
		    my $chr1=$promotercoords{$idtx1}{"chr"};
		    my $start1=$promotercoords{$idtx1}{"start"};
		    my $end1=$promotercoords{$idtx1}{"end"};

		    my $type1=$geneinfo{$gene1}{"biotype"};
		    
		    foreach my $idtx2 (keys %{$overlap{$idother}}){
			my $gene2=$txcoords{$idtx2}{"gene"};
			my $strand2=$txcoords{$idtx2}{"strand"};
			my $chr2=$promotercoords{$idtx2}{"chr"};
			my $start2=$promotercoords{$idtx2}{"start"};
			my $end2=$promotercoords{$idtx2}{"end"};

			my $type2=$geneinfo{$gene2}{"biotype"};
			
			if($gene1 ne $gene2){
			    print $output $idtx1."\t".$gene1."\t".$type1."\t".$chr1."\t".$start1."\t".$end1."\t".$strand1."\t".$idtx2."\t".$gene2."\t".$type2."\t".$chr2."\t".$start2."\t".$end2."\t".$strand2."\n";

			}
		    }
		}
	    }
	}
    }
}

close($output);

print "Done.\n";

#####################################################################################
