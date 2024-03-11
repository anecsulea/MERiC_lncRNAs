use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGeneInfo{
    my $pathin=$_[0];
    my $okchromo=$_[1];
    my $geneinfo=$_[2];

    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[$header{"Gene stable ID"}];
	my $biotype=$s[$header{"Gene type"}];
	my $chr=$s[$header{"Chromosome/scaffold name"}];

	if(exists $okchromo->{$chr}){
	    $geneinfo->{$geneid}={"biotype"=>$biotype};
	}
       	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readGeneNames{
    my $pathin=$_[0];
    my $geneinfo=$_[1];
    my $genenames=$_[2];
    my $reversenames=$_[3];

    open(my $input, $pathin);
    
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[$header{"Gene stable ID"}];
	my $genename=$s[$header{"Gene name"}]; 

	if($genename ne ""){
	    if(exists $geneinfo->{$geneid}){
		if(exists $genenames->{$genename}){
		    $genenames->{$genename}{$geneid}=1;
		} else{
		    $genenames->{$genename}={$geneid=>1};
		}

		if(exists $reversenames->{$geneid}){
		    $reversenames->{$geneid}{$genename}=1;
		} else{
		    $reversenames->{$geneid}={$genename=>1};
		}
	    }
	}
       	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readSynonyms{
    my $pathin=$_[0];
    my $geneinfo=$_[1];
    my $synonyms=$_[2];
    my $reversenames=$_[3];
    
    open(my $input, $pathin);

    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $geneid=$s[$header{"Gene stable ID"}];

	if(exists $geneinfo->{$geneid}){
	    my $syn=$s[$header{"Gene Synonym"}];

	    if($syn ne ""){
		if(exists $synonyms->{$syn}){
		    $synonyms->{$syn}{$geneid}=1;
		} else{
		    $synonyms->{$syn}={$geneid=>1};
		}
		
		if(exists $reversenames->{$geneid}){
		    $reversenames->{$geneid}{$syn}=1;
		} else{
		    $reversenames->{$geneid}={$syn=>1};
		}
	    }
	} 
	
	$line=<$input>;
    }

    close($input);

    ## remove ambiguous synonyms

    my %ambiguous;

    my @allsyn=keys %{$synonyms};

    foreach my $syn (@allsyn){
	my $nbg=keys %{$synonyms->{$syn}};

	if($nbg>=2){
	    $ambiguous{$syn}=1;
	    delete $synonyms->{$syn};
	}
    }

    my $nbambiguous=keys %ambiguous;

    print "There are ".$nbambiguous." ambiguous synonyms, we removed them.\n";
}

##############################################################

sub readForbiddenGenes{
    my $pathfg=$_[0];
    my $genelist=$_[1];

    open(my $input, $pathfg);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split(" ", $line);
	my $gene=$s[0];

	if($gene ne ""){
	    $genelist->{$gene}=1;
	}
	
	$line=<$input>;
    }
    
    close($input);

    my $nbf=keys %{$genelist};

    print "There are ".$nbf." forbidden genes.\n";
}

##############################################################

sub parsePubMedResults{
    my $pathin=$_[0];
    my $results=$_[1];

    my @lines;
    open(my $input, $pathin);
    my $line=<$input>;

    my $currentpmid="NA";
    my $lastfield="NA";
    
    while($line){
	chomp $line;

	my $prefix=substr $line,0,5;
	my $restofline=substr $line, 6;

	# if($prefix eq "     "){
	#     print "found empty prefix, last prefix ".$lastfield."\n";
	# }

	## keep track of the last info we saw

	    
	if($prefix eq "PMID-"){
	    $currentpmid=$restofline;

	    ## initialize PubMed id
	    $results->{$currentpmid}={"abstract"=>"", "title"=>"", "journal"=>"",  "type"=>"", "year"=>""};
	}

	## abstract
	if($prefix eq "AB  -"){
	    $results->{$currentpmid}{"abstract"}=$restofline;
	}

	## abstract can be on multiple lines
	if($prefix eq "     " && $lastfield eq "AB  -"){
	    $results->{$currentpmid}{"abstract"}.=$restofline;
	}

	## title
	if($prefix eq "TI  -"){
	    $results->{$currentpmid}{"title"}=$restofline;
	}

	## title can be on multiple lines
	if($prefix eq "     " && $lastfield eq "TI  -"){
	    $results->{$currentpmid}{"title"}.=$restofline;
	}
	
	## entry type
	if($prefix eq "PT  -"){
	    $results->{$currentpmid}{"type"}=$restofline;
	}

	## journal name
	if($prefix eq "JT  -"){
	    $results->{$currentpmid}{"journal"}=$restofline;
	}

	## year of publication
	if($prefix eq "DP  -"){
	    my @s=split(" ", $restofline);
	    my $year=$s[0];
	    $results->{$currentpmid}{"year"}=$year;
	}
	
	if($prefix ne "     "){
	    $lastfield=$prefix;
	}
	
	$line=<$input>;
    }
    
    close($input);
    
}

##############################################################

sub extractAssociatedGenes{
    my $results=$_[0];
    my $genenames=$_[1];
    my $synonyms=$_[2];
    my $forbidden=$_[3];
    my $nboccurrences=$_[4];

    foreach my $pmid (keys %{$results}){

	$results->{$pmid}{"geneids"}={};

	my $names="";
	my $ids="";
	
	my $abstract=$results->{$pmid}{"abstract"};
	my @words=split(" ", $abstract);

	my %genes;
	
	
	foreach my $word (@words){
	    ## remove parantheses
	    $word =~ s/\(//;
	    $word =~ s/\)//;
	    ## square brackets
	    $word =~ s/\[//;
	    $word =~ s/\]//;
	    ## remove commas
	    $word =~ s/\,//;
	    $word =~ s/\;//;

	    # $word=uc $word;

	    if(!exists $forbidden->{$word}){
		if(exists $genenames->{$word} || exists $synonyms->{$word}){
		    $genes{$word}=1;
		}
	    }
	}

	foreach my $g (keys %genes){
	    $names.=$g.";";

	    if(exists $nboccurrences->{$g}){
		$nboccurrences->{$g}++;
	    } else{
		$nboccurrences->{$g}=1;
	    }

	    if(exists $genenames->{$g}){
		foreach my $geneid (keys %{$genenames->{$g}}){
		    $results->{$pmid}{"geneids"}{$geneid}=1;
		}
	    } else{
		if(exists $synonyms->{$g}){
		    foreach my $geneid (keys %{$synonyms->{$g}}){
			$results->{$pmid}{"geneids"}{$geneid}=1;
		    }
		} 
	    }
	}
    }
}

##############################################################

sub writeOutput{
    my $results=$_[0];
    my $reversenames=$_[1];
    my $pathout=$_[2];

    open(my $output, ">".$pathout);
    
    foreach my $pmid (keys %{$results}){
	print $output "PMID: ".$pmid."\n";
	print $output "Title: ".$results->{$pmid}{"title"}."\n";
	print $output "Journal: ".$results->{$pmid}{"journal"}."\n";
	print $output "Year: ".$results->{$pmid}{"year"}."\n";
	print $output "Abstract: ". $results->{$pmid}{"abstract"}."\n";

	my @keywords=("long noncoding", "long non-coding", "Long non-coding", "Long noncoding", "lncRNA", "lincRNA", "lncRNAs", "lincRNAs");

	my $haslnc="False";

	foreach my $kw (@keywords){
	    my $hl=grep(/$kw/, $results->{$pmid}{"abstract"});

	    if($hl>0){
		$haslnc="True";
		last;
	    }
	}
	
	my $gstr="Associated genes: ";

	foreach my $gene (keys %{$results->{$pmid}{"geneids"}}){
	    $gstr.=$gene. " (";

	    foreach my $name (keys %{$reversenames->{$gene}}){
		$gstr.=$name.",";
	    }

	    chop $gstr;
	    $gstr.="); ";
	}
	
	print $output $gstr."\n";

	print $output "Mentions lncRNAs: ".$haslnc."\n";
	print $output "\n";
    }

    close($output);
}


##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts article information from list of PubMed results.\n";
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

$parameters{"pathPubMedResults"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathGeneNames"}="NA";
$parameters{"pathSynonyms"}="NA";
$parameters{"pathForbiddenGenes"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathPubMedResults", "pathGeneInfo", "pathGeneNames", "pathSynonyms", "pathForbiddenGenes", "pathOutput");

my %numericpars;

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

print "Parsing PubMed results...\n";

my %results;

my $path=$parameters{"pathPubMedResults"};

parsePubMedResults($path, \%results);

my $nbentries=keys %results;

print "Found ".$nbentries." entries.\n";

print "Done.\n";

#####################################################################

print "Reading gene information...\n";

my %okchromo;

for(my $i=1; $i<=22; $i++){
    $okchromo{$i.""}=1;
}

$okchromo{"X"}=1;
$okchromo{"Y"}=1;
$okchromo{"MT"}=1;

my $nbok=keys %okchromo;

print "We accept genes on ".$nbok." chromosomes.\n";

my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%okchromo, \%geneinfo);

my $nbgenes=keys %geneinfo;

print "There are ".$nbgenes." genes.\n";

print "Done.\n";
    
#####################################################################

print "Reading gene names...\n";

my %genenames;
my %reversenames;

readGeneNames($parameters{"pathGeneNames"}, \%geneinfo, \%genenames, \%reversenames);

my $nbnames=keys %genenames;

print "Found ".$nbnames." different gene names.\n";

print "Done.\n";

#####################################################################

print "Reading gene name synonyms...\n";

my %synonyms;

readSynonyms($parameters{"pathSynonyms"}, \%geneinfo, \%synonyms, \%reversenames);

my $nbsyn=keys %synonyms;

print "Found ".$nbsyn." synonyms.\n";

print "Done.\n";

#####################################################################

print "Reading forbidden genes...\n";

my %forbidden;

readForbiddenGenes($parameters{"pathForbiddenGenes"}, \%forbidden);

print "Done.\n";

#####################################################################

print "Extracting gene names associated with articles...\n";

my %nbocc;

extractAssociatedGenes(\%results, \%genenames, \%synonyms, \%forbidden, \%nbocc);

# statistics for numbers of genes
# my %nbgenesocc;

# foreach my $gene (keys %nbocc){
#     my $nb=$nbocc{$gene};

#     if(exists $nbgenesocc{$nb}){
# 	push(@{$nbgenesocc{$nb}}, $gene);
#     } else{
# 	$nbgenesocc{$nb}=[$gene];
#     }
# }

# my @nbocc=sort {$a<=>$b} (keys %nbgenesocc);

# foreach my $nb (reverse @nbocc){
#     print $nb." ".join(", ", @{$nbgenesocc{$nb}})."\n";
# }

print "Done.\n";

#####################################################################
   
print "Writing output...\n";

writeOutput(\%results, \%reversenames, $parameters{"pathOutput"});

print "Done.\n";

#####################################################################
#####################################################################
