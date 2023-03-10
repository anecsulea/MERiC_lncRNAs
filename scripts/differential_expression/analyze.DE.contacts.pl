use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###################################################################
###################################################################

sub readDifferentiallyExpressedGenes{
    my $pathin=$_[0];
    my $maxFDR=$_[1];
    my $minFC=$_[2];
    my $signifgenes=$_[3];

    my $thresholdfc=log($minFC)/log(2);

    open(my $input, $pathin);

    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    $header{"GeneID"}=0; ## formatted by R with row names
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i+1;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $geneid=$s[$header{"GeneID"}];
	my $fdr=$s[$header{"padj"}]+0.0;
	my $logfc=$s[$header{"log2FoldChange"}]+0.0;

	if($fdr<$maxfdr && abs($logfc)>=$thresholdfc){
	    $signifgenes->{$geneid}=1;
	}
	
	$line=<$input>;
    }

    close($input);
}

###################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script analyzes DE genes with respect to PCHi-C contacts. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###################################################################
###################################################################

my %parameters;

$parameters{"pathDifferentialExpression"}="NA";
$parameters{"maxFDR"}="NA";
$parameters{"minFC"}="NA";
$parameters{"pathPCHiCInteractions"}="NA";
$parameters{"minDistance"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathDifferentialExpression", "maxFDR", "minFC", "pathPCHiCInteractions", "minDistance", "pathOutput");

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

print "Reading DE genes...\n";

my $minfc=$parameters{"minFC"}+0.0;
print "minimum fold change: ".$minfc."\n";

my $maxfdr=$parameters{"maxFDR"}+0.0;
print "max FDR: ".$maxfdr."\n";

my %signifgenes;
readDifferentiallyExpressedGenes($parameters{"pathDifferentialExpression"}, $maxfdr, $minfc, \%signifgenes);

my $nbsignif=keys %signifgenes;

print "There are ".$nbsignif." significantly expressed genes.\n";

print "Done.\n";

#####################################################################




#####################################################################
