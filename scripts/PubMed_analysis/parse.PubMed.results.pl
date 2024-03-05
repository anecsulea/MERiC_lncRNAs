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
    my $type=$_[1];
    my $results=$_[2];

    my $nbeabstract=4;
    my $nbetitle=1;
    my $indexjournal=0;
    my $splitjournal=1;
    
    if($type eq "retraction"){
	$nbeabstract=6;
	$nbetitle=2;
	$indexjournal=2;
	$splitjournal=0;
    }

    print "Abstract should be the ".$nbeabstract."th entry in this file\n";
    print "Title should be on line ".$nbetitle." in this file\n";

    
    my @lines;
    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	push(@lines, $line);
	$line=<$input>;
    }
    
    close($input);
    
    my %entryindexes;
    my %duplicated;

    my $currententry=1;

    for(my $i=0; $i<@lines; $i++){
	my $line=$lines[$i];
	
	my @s=split("\\. ", $line);
	my $prefix=$s[0]+0;

	if($prefix==$currententry){
	    $entryindexes{$currententry}=$i;
	    $currententry++;
	}
    }

    my $nbentries=$currententry-1;

    for(my $entry=1; $entry<=$nbentries; $entry++){
	my $firstindex=$entryindexes{$entry};
	my $lastindex="NA";
	
	if($entry==$nbentries){
	    $lastindex=@lines-1;
	} else{
	    $lastindex=$entryindexes{$entry+1}-1;
	}

	my $firstline=$lines[$firstindex+$indexjournal];
	my @s=split("\\.", $firstline);
	my $journal=$s[$splitjournal];
	my @date=split(" ", $s[$splitjournal+1]);
	my @yearinfo=split(";", $date[0]);
	my @yyinfo=split(":", $yearinfo[0]);
	my @yyyinfo=split("-", $yyinfo[0]);
	my $year=$yyyinfo[0];

	my $abstract="";
	my $title="";
	my $pmid="NA";
	
	## count number of empty lines in entry

	my $nbemptylines=0;
	
	for(my $i=$firstindex; $i<=$lastindex; $i++){
	    if($lines[$i] eq ""){
		$nbemptylines++;
	    }
	}

	if($nbemptylines>=6){
	    my $nbabstractlines=0;
	   	    
	    my $currentemptyline=0;
	    
	    for(my $i=$firstindex; $i<=$lastindex; $i++){
		
		if($lines[$i] eq ""){
		    $currentemptyline++;
		} else{
		    if($currentemptyline==$nbetitle){
			if($title eq ""){
			    $title=$lines[$i];
			} else{
			    $title.=" ".$lines[$i];
			}
		    }
		    my @st=split(":", $lines[$i]);
		    my $prefix=$st[0];

		    if($currentemptyline==$nbeabstract){
			if($prefix ne "DOI" && $prefix ne "PMCID" && $prefix ne "PMID"){
			    $nbabstractlines++;

			    if($abstract eq ""){
				$abstract=$lines[$i];
			    } else{
				$abstract.=" ".$lines[$i];
			    }
			}
		    }

		    if($prefix eq "PMID"){
			my @tt=split(" ", $st[1]);
			$pmid=$tt[0];
		    }
		}
	    }

	    if($nbabstractlines>2 && $pmid ne "NA"){
		my @st=split(":", $abstract);
		my $prefix=$st[0];

		if($prefix ne "Author information"){
		    if(exists  $results->{$pmid}){
			if($results->{$pmid}{"title"} ne $title || $results->{$pmid}{"journal"} ne $journal || $results->{$pmid}{"year"} ne $year || $results->{$pmid}{"abstract"} ne $abstract){ 
			    print "Weird ! we already saw this PubMed id: ".$pmid.", with different info\n";
			    print "year ".$results->{$pmid}{"year"}." ".$year."\n";
			    print "abstract ".$results->{$pmid}{"abstract"}." ".$abstract."\n";
			    print "title ".$results->{$pmid}{"title"}." ".$title."\n";
			    print "journal ".$results->{$pmid}{"journal"}." ".$journal."\n";
			    $duplicated{$pmid}=1;
			}
		    } else{
			$results->{$pmid}={"title"=>$title, "journal"=>$journal, "year"=>$year, "abstract"=>$abstract};
		    }
		}
	    }
	}
    }

    my $nbres=keys %{$results};

    print "There are ".$nbres." different PMIDs.\n";
    
    my $nbdupli=keys %duplicated;

    foreach my $id (keys %duplicated){
	delete $results->{$id};
    }
    
    print "There were ".$nbdupli." duplicated entries.\n";
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

$parameters{"citationType"}="NA";
$parameters{"pathPubMedResults"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"pathGeneNames"}="NA";
$parameters{"pathSynonyms"}="NA";
$parameters{"pathForbiddenGenes"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("citationType","pathPubMedResults", "pathGeneInfo", "pathGeneNames", "pathSynonyms", "pathForbiddenGenes", "pathOutput");

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

my @paths=split(",", $parameters{"pathPubMedResults"});

my $type=$parameters{"citationType"};

foreach my $path (@paths){
    print "Reading info from ".$path."\n";
    parsePubMedResults($path, $type, \%results);
}

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
