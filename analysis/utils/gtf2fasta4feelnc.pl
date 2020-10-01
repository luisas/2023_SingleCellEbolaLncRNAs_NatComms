use lib "/home/luisas/Desktop/FEELnc-master/lib";



#!/usr/bin/perl

#
# Modification by V.Wucher april 16 2015:
# 	Modification of the predicting method: use now random forest (package R randomForest)
#

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use List::Util qw( min max );

# lib directory: {FEELnc github directory}/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;
use RandomForest;
use ExtractCdnaOrf;
use RNAshuffle;

# Program name
my $progname = basename($0);

# Variables
my $mRNAfile = "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_known_proteincoding_test.gtf";
my $genome = "/home/luisas/Desktop/cluster/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa";
my $lncRNAfile = undef;
my %biotype;
my $man        = 0;
my $help       = 0;
my $verbosity  = 1;
# my $outputlog;
my $numtx    = undef; # number of mRNAs and lncRNAs tx for training separate by a ','. undef or '.'  for all transcripts
my $minnumtx = 100;   # Min number of tx for training (a too small value will result in a bad learning)


# VW Add a variable to get the kmer size which are used to get the kmer scores
my $kmerList = '1,2,3,6,9,12';

# VW Add a variable to keep tmp file, default don't keep
my $keepTmp = 0;

# VW If random forest (rf/RF) cutoff is defined, no need to compute it on TP lncRNA and mRNA
#    and a $speThres for the thresolds on mRNA and lncRNA specificity
my $rfcut        = undef;
my $speThres     = undef;
my @speThresList = undef;

# VW Add option to select the calculate orf for learning and test data sets
my $orfTypeLearn = 3;
my $orfTypeTest  = 3;

# VW Add an option to specify the output directory, default current directoryand an out name
my $outDir  = "./feelnc_codpot_out";
my $outName = "";

# VW Add an option to change the number of trees use in random forest
my $nTree = 500;

# VW Add an option to fixe the seed
my $seed = 1234;

# VW Add nbr proc
my $proc = 1;

# VW Add a percentage to get two learning file
my $perc = 0.1;

# VW Add the mode option to known if FEELnc need to extract intergenic sequences or shuffle mRNAs sequences
# if there is no input lncRNAs
my $mode = "";

# Intergenic extraction:
my $maxTries   = 10;
my $maxN       = 5;
my $sizecorrec = 0.75; # a float value between 0 and 1

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
    'a|mRNAfile=s'   => \$mRNAfile,
    'g|genome=s'     => \$genome,
    'n|numtx=s'      => \$numtx,
    'b|biotype=s'    => \%biotype,
    'r|rfcut=f'      => \$rfcut,
    'spethres=s'     => \$speThres,
    'k|kmer=s'       => \$kmerList,
    'm|mode=s'       => \$mode,
    's|sizeinter=f'  => \$sizecorrec,
    'learnorftype=i' => \$orfTypeLearn,
    'testorftype=i'  => \$orfTypeTest,
    'ntree=i'        => \$nTree,
    'outdir=s'       => \$outDir,
    'o|outname=s'    => \$outName,
    'percentage=f'   => \$perc,
    'keeptmp'        => \$keepTmp,
    'v|verbosity=i'  => \$verbosity,
    'p|processor=i'  => \$proc,
    'seed=i'         => \$seed,
    'help|?'         => \$help,
    'man'            => \$man
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("- Error: Cannot read your input annotation file '$mRNAfile'...\nFor help, see:\n$progname --help\n") unless( -r $mRNAfile);
if (defined $rfcut){
    pod2usage ("- Error: --rfcut option '$rfcut' should be a float between 0 and 1 [0-1] \n") unless ($rfcut >= 0 and $rfcut <= 1);
}
pod2usage ("- Error: --sizecorrec option (ratio between mRNAs sequence lenghts and intergenic non coding sequence lenghts) '$sizecorrec' should be a float between 0 and 1 [0-1]\n") unless ($sizecorrec >= 0 and $sizecorrec <= 1);
pod2usage ("- Error: --orfTypeLearn option '$orfTypeLearn' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information)\n") unless ($orfTypeLearn==0 || $orfTypeLearn==1 || $orfTypeLearn==2 || $orfTypeLearn==3 || $orfTypeLearn==4);
pod2usage ("- Error: --orfTypeTest option '$orfTypeTest' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information)\n") unless ($orfTypeTest==0 || $orfTypeTest==1 || $orfTypeTest==2 || $orfTypeTest==3 || $orfTypeTest==4);
# pod2usage ("- Error: --outDir option '$outDir' is not a directory or it does not exist \n") unless (-d $outDir);
pod2usage ("- Error: --nTree option '$nTree' should be strictly positive\n") unless ($nTree > 0);
pod2usage ("- Error: --rfcut and --spethres specified, only one of the two options can be used (default one threshold defined on a 10-fold cross-validation)\n") if((defined $rfcut) && (defined $speThres));
pod2usage ("- Error: -p/--processor option '$proc' should be a positive integer\n") unless ($proc >= 1);
pod2usage ("- Error: --percentage option '$perc' should be a number in ]0;1[\n") unless ($perc>0 && $perc<1);

# If -m|--mode is set, check if the value is either intergenic or shuffle
pod2usage ("- Error: the value set for the -m|--mode option is not valide. It need to be either 'shuffle' or 'intergenic'\n") if ( ($mode ne "") && ($mode !~ /shuffle|intergenic/) );

# If no lncRNA sequences and any option between --intergenic and --shuffle have been choosen, quit




# Check the max kmersize
my @kmerTable = split(/,/,$kmerList);
my $kmerMax   = max @kmerTable;
pod2usage ("- Error: \$kmerList option '$kmerList' is not valid. One of the size is stricly greater than '15' (see --help for more details)\n") unless ($kmerMax <= 15);

# Check threshold values for the mRNAs and lncRNAs specificity
if((defined $speThres))
{
    @speThresList = split(/,/, $speThres);
    if(@speThresList!=2)
    {
	pod2usage ("- Error: --speThres option '$speThres' should be a list of two value separated by a ',' (see --help for more details)\n");
    }
    if(($speThresList[0]<=0) || ($speThresList[0]>=1) || ($speThresList[1]<=0) || ($speThresList[1]>=1))
    {
	pod2usage ("- Error: one value of --speThres option '$speThres' is equal or greater than 1 or equal or lesser than 0, should be in ]0,1[ (see --help for more details)\n");
    }
}



# Check the number of mRNAs ($numtxCod) and lncRNAs ($numtxNon) transcripts use for learning
my $numtxCod = undef;
my $numtxNon = undef;

if(defined $numtx)
{
    if(($numtx =~ tr/,//) != 1)
    {
	pod2usage ("- Error: --numtx option '$numtx' should be a list of two value separated by a ',' (see --help for more details)\n");
    }
    else
    {
	($numtxCod,$numtxNon) = split(/,/, $numtx);
    }

    if($numtxCod =~ /^\d+$/ &&  $numtxCod < $minnumtx)
    {
	print $numtxCod."     ".$minnumtx."\n";
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif(($numtxCod !~ /^\d+$/) && $numtxCod ne "")
    {
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif($numtxCod eq "")
    {
	$numtxCod = undef;
    }

    if($numtxNon =~ /^\d+$/ &&  $numtxNon < $minnumtx)
    {
	print $numtxNon."     ".$minnumtx."\n";
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif(($numtxNon !~ /^\d+$/) && $numtxNon ne "")
    {
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif($numtxNon eq "")
    {
	$numtxNon = undef;
    }
}




# Create the directory for temporary files with the job id ($$) and the temporary name
my $outTmp = "/home/luisas/Desktop";
my $outName=basename($mRNAfile);
my $nameTmp = $outTmp."/".$outName;


# Initialize the seed
srand($seed);

# If $numtx is undef, then learning on all transcripts, can be long so print a warning...
print STDERR "You do not have specified a maximum number mRNAs transcripts for the training. Use all the annotation, can be long...\n"  if(!defined $numtxCod);
print STDERR "You do not have specified a maximum number lncRNA transcripts for the training. Use all the annotation, can be long...\n" if(!defined $numtxNon);



# Define fasta file names
my $codFile    = $nameTmp.".coding_rna.fa";
my $codOrfFile = $nameTmp.".coding_orf.fa";
my $nonFile    = $nameTmp.".noncoding_rna.fa";
my $nonOrfFile = $nameTmp.".noncoding_orf.fa";
my $testFile    = $nameTmp.".test_rna.fa";
my $testOrfFile = $nameTmp.".test_orf.fa";



##########################################################
# mRNA file
#######
# add a refhash that will contain the mRNA ID that passed cDNA and ORF steps
# Will be used to checked for randomization
my $ref_cDNA_passed;
my $refmrna;


($ref_cDNA_passed, $refmrna) = ExtractCdnaOrf::CreateORFcDNAFromGTF($mRNAfile, $codFile, $codOrfFile, $numtxCod, $minnumtx, $genome, 'exon,CDS,stop_codon,start_codon', \%biotype, $orfTypeLearn, $verbosity, $kmerMax);
