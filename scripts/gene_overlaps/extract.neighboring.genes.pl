use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readBlocks{
    my $pathin=$_[0];
    my $blocks=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $gene=$s[0];
	my $chr=$s[2];
	my $start=$s[3];
	my $end=$s[4];
	my $strand=$s[5];

	if(!(exists $blocks->{$gene})){
	    $blocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[], "end"=>[]};
	}

	push(@{$blocks->{$gene}{"start"}}, $start);
	push(@{$blocks->{$gene}{"end"}}, $end);
	
	$line=<$input>;
    }

    close($input);
}

###################################################################

sub extractGeneCoords{
    my $blocks=$_[0];
    my $genecoords=$_[1];

    foreach my $gene (keys %{$blocks}){
	my $chr=$blocks->{$gene}{"chr"};
	my $strand=$blocks->{$gene}{"strand"};

	my $start=min @{$blocks->{$gene}{"start"}};
	my $end=max @{$blocks->{$gene}{"end"}};

	$genecoords->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>$start, "end"=>$end};
    }
}

###########################################################################

sub orderGeneCoords{
    my $coords=$_[0];
    my $ordered=$_[1];
    
    my %hashpos;

    foreach my $gene (keys %{$coords}){
	my $chr=$coords->{$gene}{"chr"};
	my $strand=$coords->{$gene}{"strand"};
	my $start=$coords->{$gene}{"start"};
	my $end=$coords->{$gene}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"gene"}},$gene);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		push(@{$hashpos{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$start}={"gene"=>[$gene],"end"=>[$end],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"gene"=>[$gene],"end"=>[$end],"strand"=>[$strand]}};
	}
    }
   

    foreach my $chr (keys %hashpos){
	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];
		my $gene=${$hashpos{$chr}{$start}{"gene"}}[$i];
		
		if(!(exists $ordered->{$chr})){
		    $ordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};
		}

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"strand"}}, $strand);	
		push(@{$ordered->{$chr}{"id"}}, $gene);
	    }
	}
    }
}

######################################################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $type=$_[3];
    my $overlap=$_[4];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"start"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$coords1->{$chr}{"start"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"end"}}[$i]+$margin;

		my $strand1=${$coords1->{$chr}{"strand"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];
			
		my $j=$firstj;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $strand2=${$coords2->{$chr}{"strand"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];
		    
		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense") || ($type eq "any")){
			
			my $M=max($start1,$start2);
			my $m=min($end1,$end2);
			
			if($M<=$m){
			    my $fr1=($m-$M+1)/($end1-$start1+1);
			    my $fr2=($m-$M+1)/($end2-$start2+1);
			    
			    if(exists $overlap->{$id1}){
				$overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2};
			    }
			    else{
				$overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2}};
			    }
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}

###################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts neighboring genes. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###################################################################

my %parameters;

$parameters{"pathExonBlocks"}="NA";
$parameters{"windowSize"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonBlocks", "windowSize", "pathOutput");

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

print "Reading exon block coordinates...\n";

my %exonblocks;

readBlocks($parameters{"pathExonBlocks"}, \%exonblocks);

my $nbg=keys %exonblocks;

print "There are ".$nbg." genes in total.\n";

print "Done.\n";

#####################################################################

print "Extracting gene coordinates...\n";

my %genecoords;

extractGeneCoords(\%exonblocks, \%genecoords);  

print "Done.\n";

#####################################################################

print "Ordering coordinates...\n";

my %orderedgenes;

orderGeneCoords(\%genecoords, \%orderedgenes); 

print "Done.\n";

#####################################################################

print "Extracting neighboring genes...\n";

my $window=$parameters{"windowSize"}+0;

print "Window size: ".$window."\n";

my %ovgenessense;
extractOverlap(\%orderedgenes, \%orderedgenes, $window, "sense", \%ovgenessense);

my $nbov=keys %ovgenessense;
print "There are ".$nbov." genes with sense genic overlap.\n";


my %ovgenesantisense;
extractOverlap(\%orderedgenes, \%orderedgenes, $window, "antisense", \%ovgenesantisense);

my $nbov=keys %ovgenesantisense;
print "There are ".$nbov." genes with antisense genic overlap.\n";

print "Done.\n";

#####################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tNeighborID\tType\n";

foreach my $gene (keys %genecoords){
   
    ## we analyze overlap between whole genes
    
    if(exists $ovgenessense{$gene}){
	foreach my $og (keys %{$ovgenessense{$gene}}){
	    if($og ne $gene){
		print $output $gene."\t".$og."\tsense\n";
	    }
	}
    }

    if(exists $ovgenesantisense{$gene}){
	foreach my $og (keys %{$ovgenesantisense{$gene}}){
	    if($og ne $gene){
		print $output $gene."\t".$og."\tantisense\n";
	    }
	}
    }
}


print "Done.\n";

#####################################################################
