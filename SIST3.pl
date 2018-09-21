####################
# SIST: Separate Accugenomics Spike-In reads (3 bins version)
# Version 3.01
# 
# #### Updates: ####
# 09/14/2018: fix some bugs and adding error message of missing picard Ver. 3.01)
# 09/12/2018: Major Change: (three bins stratey) adding one category called "suspected reads" which includes low quality reads. (forked version: Ver. 3.00)
#			  the version jumps to 3.00 in order to the three bins strategy. the two version will update separately.
# 08/16/2018: using git-hub for version control (update information depreciated)
# 08/02/2018: provide an option to use customized spike-in reference. (Ver. 1.10 -> 1.20)
# 06/09/2018: Support BAM output (to include additional read information)
# 01/04/2018: Support BAM input (Ver. 1.01 -> 1.10)
# 12/02/2017: Fix a Bug about read head detection in FASTQ file generation step. (ver. 1.00 -> 1.01)
# #### End of Updates ###
#
# Contact leihong.wu@fda.hhs.gov for further assistance.
#
####################

use Getopt::Long;
use strict;

########################################################
# USAGE
my $USAGE =<<USAGE;
    Usage:
      SIST3.pl [-fastq_1=<fastq.gz> -fastq_2=<fastq.gz>] [-bam=<bam>] -O=<output_prefix> [options]
    
    Required inputs:
      -fastq_1  paired-end reads 1 (required if no bam input)
      -fastq_2  paired-end reads 2 (required if paired-end reads)
    
      -bam  input file as bam format (required if no fastq input)
    
      -O    Output file header 
    
    Common options:
      -ref [prefix]: prefix of spike-in reference, (defalt: "Refs/Accugenomics_Spikein_V2");
	  -type [match, all, gzip] 
        match: only extract the spike-in reads 
        all:   generate fastq file for spike-in reads and remain origin reads (default)
        gzip:  generate gz file instead of fastq file
      -keepBam : when the input is bam file and keepBam is specified, output would keep the sam format and not change to fastq; 
      -umi: flag for umi reads;
      
    Other options:
      -t [num]: threads used in BWA mem (default: 8);
      -bwa_b [num]: mismatch penalty used in BWA mem (default: 6);
      -bwa_o [num,num]: indel penalty used in BWA mem (default: [20,20]);
      -min_length [num]: minimum match length for a read (default: 70);
      -r_hg : human reference genome (default: Refs/genome.fa);
      -picard : picard path (default: Refs/picard.jar);
      
    -help:  Prints out this helpful message
    
USAGE
#
######################################################

my $SAMPLE_1 = '';
my $SAMPLE_2 = '';
my $BAM = '';
my $hq_header = '';

my $threads = 4;
my $type = 'all';
my $BWA_B= 6; # penalty for mismatch. default = 4 in BWA MEM for short reads
my $BWA_O= '20,20';  # penalty for indels. default = 6, high value because indels are highly not expected in spike-in sequence.
my $map_score = 70; # at least mapped 70 bps in consecutive to the spike-in Reference. 
my $single_mode = 'false';
my $samtools_path = 'samtools';
my $test_mode = 'false';
my $ref_prefix = 'Refs/ISREF';
my $picard_file = 'Refs/picard.jar';
my $mismatch = 1; # what level of mismatch are allowed.

#### TESTING MODE ONLY #####
# my $samtools_path = '/storage2/lwu/SEQC2/samtools-1.8/build/bin/samtools'; # specific samtools path for test use.
############################

GetOptions (   
        "fastq_1=s" => \$SAMPLE_1,
        "fastq1=s"      => \$SAMPLE_1,
        "fq_1=s"        => \$SAMPLE_1,
        "fq1=s"     => \$SAMPLE_1,
        "f1=s"      => \$SAMPLE_1,
        "1=s"       => \$SAMPLE_1,
        
        "fastq2=s"      => \$SAMPLE_2,
        "fastq_2=s"     => \$SAMPLE_2,
        "fq_2=s"    => \$SAMPLE_2,
        "fq2=s"     => \$SAMPLE_2,
        "f2=s"      => \$SAMPLE_2,
        "2=s"       => \$SAMPLE_2,
             
        "bam=s"     => \$BAM,
             
        "O=s"       => \$hq_header,
        "o=s"       => \$hq_header,
             
        "ref=s"     => \$ref_prefix,
        "r=s"       => \$ref_prefix,
        "R=s"       => \$ref_prefix,
        "Ref=s"     => \$ref_prefix,
		
		"maxmismatch=i" => \$mismatch,
		"m=i" => \$mismatch,
		
		"picard=s"	=> \$picard_file,
             
        "type=s"    => \$type, 
        "t=s"       => \$threads,
        "bwa_b=s"   => \$BWA_B,
        "bwa_o=s"   => \$BWA_O,
        "min_length=s"  => \$map_score,
        q(help)     => \my $help,
        q(h)        => \my $help,
        q(keepBam)  => \my $keepBam,
		q(umi)  	=> \my $umi_tag,
        # ONLY USED FOR DEV PURPOSES;
        "test_mode=s"   => \$test_mode
);            

## TEST_MODE                     
($SAMPLE_1, $SAMPLE_2, $BAM, $hq_header) = get_test_data($test_mode) if ($test_mode ne 'false');

## Code Section 1. Pre-Check Begin ##
# c1.1 Help message
if ($help || (!$SAMPLE_1 && !$BAM)) {
    # print "\tSomething Wrong\n";
    print "$USAGE\n";
    exit 0;
}

# c1.2 Check analysis type (-type)
if ($type ne 'match' && $type ne 'all' && $type ne 'gzip' ){
    print('Wrong type! should be one in [match, all, gzip] '."\n");
    exit 0;
}

# c1.3 Check Picard (only works when doing insert size measurement)
#if (! -f 'Refs/picard.jar' && $eval){
#   print('Missing Picard. create a link of picard  Refs/picard.jar'."\n");
#   exit 0;
#}

# c1.5 Check Input Files Availability (-f1, -f2, -bam)
if ( (!$SAMPLE_2)&& !$BAM){
    $single_mode = 'true'; # single-end reads are used in this study
    print "Detect current analysis uses single-ended reads; 
        if not, please stop and double check whether fastq_2 is correct. \n";
}

if ( (!$BAM||! -f $BAM) && (!$SAMPLE_1||! -f $SAMPLE_1 )){  
    print $BAM ;
    print " Missing input files!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

if ( (!$BAM||! -f $BAM) && ($single_mode eq 'false') && (!$SAMPLE_2||! -f $SAMPLE_2 )){ 
    print " Missing paired input files!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

if ($BAM){
    if ($BAM !~ /\.bam$/ ){
        print "Not supported BAM input file format! Must ended with .bam \n";
        exit 0;
    }
}

if ($SAMPLE_1){
    if ($SAMPLE_1 !~ /\.fa*s*t*q$/ && $SAMPLE_1 !~ /\.fa*s*t*q\.gz$/){
        print "Not supported Fastq file format for Fastq_1! Must ended with .fastq[.gz] or .fq[.gz] \n";
        exit 0;
    }
}

if ($SAMPLE_2){
    if ($SAMPLE_2 !~ /\.fa*s*t*q$/ && $SAMPLE_2 !~ /\.fa*s*t*q\.gz$/){
        print "Not supported Fastq file format for Fastq_2! Must ended with .fastq[.gz] or .fq[.gz] \n";
        exit 0;
    }
}

# c1.6 Check Output Header
if (!$hq_header){   
    print " Missing output folder!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

# c1.7 Check human reference genome and picard
if (! -f $picard_file){
	print "Picard is needed for the analysis, please make sure it is in the Refs/ folder as Refs/picard.jar"."\n";
	exit 0;
}

# c1.8 Create tmp folder if not existed.
my $tmp_dir = 'tmp';
mkdir($tmp_dir) unless(-d $tmp_dir);

## Pre-Check End ##


## Code Section 2. Variable Preparation Begin ##
# c2.1 Create output files
my $hq_file_is  = $hq_header.".IS.txt";
my $hq_file_ntm  = $hq_header.".NTM.txt";
my $hq_file_ntns  = $hq_header.".NTNS.txt";
my $hq_file_ntum  = $hq_header.".NTUM.txt";
my $hq_file_sus = $hq_header.".SUS.txt";
my $sam_file_sus = $hq_header.".SUS.sam";

my $is_fastq_1  = $hq_header.".IS.1.fastq";
my $nt_fastq_1  = $hq_header.".NT.1.fastq";
my $sus_fastq_1 = $hq_header.".SUS.1.fastq";

my $is_fastq_2 = $hq_header.".IS.2.fastq";
my $nt_fastq_2 = $hq_header.".NT.2.fastq";
my $sus_fastq_2 = $hq_header.".SUS.2.fastq";

if ($keepBam){
	$is_fastq_1 = $hq_header.".IS.1.bam" ;
	$nt_fastq_1 = $hq_header.".NT.1.bam";
	$sus_fastq_1 = $hq_header.".SUS.1.bam";
}

# if output file already exists; overwrite it.
print "Warning: output files (".$hq_header.".*) already exists; will overwrite! \n" if (-f $hq_file_is);

# Accugenomics spike-in Reference

my $REF_IS = $ref_prefix.'.IS.fasta';
my $REF_NT = $ref_prefix.'.NT.fasta';
my $mut_REF =$ref_prefix.'.IS.vcf';
if (! -f $mut_REF){
    $mut_REF=prepare_ref($ref_prefix);
}

my %muts=();
open(FH,$mut_REF);
while(<FH>){
    chomp;
    my @array=split("\t");
    $muts{$array[0].":".$array[1]}=$array[2];
}
close FH;

## Variable Preparation End ##

## Code Section 3. Main Program Begin ##
###### c3.1: detect Spike-in Reads ######
## c3.1.1 initialize 
my %reads_count_total=();
my %reads_count_bwa=();

my %IS_reads =();
my %NT_reads =();
my %SUS_reads=();

my %IS_sam=();
my %NT_sam = ();
my %SUS_sam=();

# output same file
open(OFH_SAM_SUS,'>'.$sam_file_sus);

# c3.1.2 Find all mapped reads (candidate spike-in reads) from original read set.
if ($SAMPLE_1){
    if ($single_mode eq 'false'){
        # paired-end reads
        open(FH,'bwa mem -h 1 -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_IS.' '.
            $SAMPLE_1.' '.$SAMPLE_2.' 2>tmp/bwa_run.log|');
	}else{
        # single-end reads
        open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_IS.' '.
            $SAMPLE_1.' 2>tmp/bwa_run.log|');
    }
}elsif($BAM){
    # input is Bam
    open(FH,$samtools_path.' fastq '.$BAM.'|bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_IS.
            ' - 2>>tmp/bwa_run.log|');
}else{
    exit 0;
}

# read output from BWA alignment
while(<FH>){
	# print($_);
	if(/^@/){
		# print header information to sam file
		print OFH_SAM_SUS $_;     
	}else{
		my $line=$_;
		chomp($line);
		my @array=split("\t",$line);
		$reads_count_total{$array[0]}=1;

		# sequence has to be mapped and has MD code, otherwise it is considered as unmapped reads.
		# Note that, unmapped reads may include Native reads as well. 
		# In current version we don't differentiate native reads and unmapped reads.
		# They are all categorized to Native reads in the end.
		
		if ($array[2] ne '*' && $line =~ /\tMD:Z:[\d\^ACGT]+\t/){
			# Parameter initialization
			my $ins_count = 0;
			my $del_count = 0;
			
			my $left_clip=0;
			my $left_clip_hard=0;
			my $right_clip=0;
			my $right_clip_hard=0;
			my $flag_clip = 0;
			
			my $flag_length = 0;
			
			# Pre-Step : reads preparation
			$reads_count_bwa{$array[0]}=1;
			my $cigar = $array[5];
			my $sum_length = length($array[9]);

			# INDEL detection
			$ins_count = $ins_count+$& while($cigar =~ /(\d+)I/g );
			$del_count = $del_count+$& while($cigar =~ /(\d+)D/g );
			
			# FLAG for clipping : if two-sides softclips
			$left_clip = $1 if ($cigar =~ /(\d+)S\d+M/); 
			$right_clip = $1 if ($cigar =~ /\d+M(\d+)S/); 
			$left_clip_hard = $1 if ($cigar =~ /(\d+)H\d+M/); 
			$right_clip_hard = $1 if ($cigar =~ /\d+M(\d+)H/); 			
			if ($umi_tag){
				#no clipping criterion for umi reads;
				$flag_clip = 5;
			}else{
				#if not umi reads, soft-clipping criterion still works
				$flag_clip = 1 if ($right_clip >0 || $left_clip >0); # single-side softclip reads 
				$flag_clip = 2 if ($right_clip_hard >0 || $left_clip_hard >0); # single-side hardclip reads             
				$flag_clip = 5 if ($left_clip>0 && $right_clip>0); # two-side softclip (one short) reads 
				$flag_clip = 6 if ($left_clip_hard>0 && $right_clip_hard>0); # two-side hardclip (one short) reads 
				$flag_clip = 11 if ($left_clip >20 && $right_clip >20); # two-side long softclip reads 
				$flag_clip = 12 if ($left_clip_hard >20 && $right_clip_hard >20); # two-side long hardclip reads
			}
				
			## FLAG for maximum match length:
			my $mapped_seq = $sum_length-$left_clip-$right_clip-$ins_count+$del_count; 
			my $muts_thres = int($mapped_seq/50)+1;
			$flag_length = ($mapped_seq-$map_score)/10;
			
			## counting Spike-in mutation site 
			my $muts_included=0;
			foreach($array[3]..($array[3]+$mapped_seq-1)){
				my $id_tmp = $array[2].':'.$_;
				if (exists $muts{$id_tmp}){
					$muts_included ++ ;
				}
			}
					
			my $remain_spikein = 0;
			my $count_ref=0;
			my $muts_rep=$ins_count; 
			
			## find spike-in counts      
			if ($muts_included >0){ 
				# Initialize parameters
				my $md_str = $1 if $line =~ (/\tMD:Z:([\^\dACGT]+)\t/);
				my @muts_tmp = ();
				# print $md_str."\n";
				while($md_str =~ /([\^ACGT]+)/g){
					my $token = $1;
					if ($token =~ /\^/){
						push @muts_tmp, length($token)-1;
					}else{
						push @muts_tmp, length($token);
					}
					# print @muts_tmp if (length($token)>1); 
				}
				
				my @md_array = split(/[\^ACGT]+/,$md_str);
				my $pos_start = $array[3];
				pop(@md_array);
		
				# Measuring on-site mutations
				my $i = 0;
				foreach(@md_array){
					my $pos_curr = $pos_start+$_;
					$pos_start = $pos_start+$_+$muts_tmp[$i];
					my $sig=0;
					if (exists $muts{$array[2].':'.$pos_curr}){
						$count_ref++;
					}
					$i++;
				}
					
				# Measuring off-site mutations
				$remain_spikein = $muts_included-$count_ref; # this is the number of detected spike-in dinucleotides (M);
				$muts_rep = $muts_rep+scalar(@md_array)-$count_ref; 
				
				my $info = $muts_included."\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line;
				# if ($flag_clip <10 && int($array[1]) < 200){
				
				my $IS_flag = 0;
				my $NT_flag = 0;
					
				# $IS_flag = 1 if ($count_ref == 0 && $flag_length>=0);
				# $IS_flag = 1 if ($remain_spikein >= 4);
				$IS_flag = 1 if ($remain_spikein >= $muts_included-$mismatch && $muts_included > 2+$mismatch);
				$IS_flag = 1 if ($remain_spikein == $muts_included && $muts_included >=2 && $muts_included <= 2+$mismatch);
				$NT_flag = 1 if ($count_ref >= $muts_included-$mismatch  && $muts_included > 2+$mismatch );
				$NT_flag = 1 if ($count_ref == $muts_included  && $muts_included >=2 && $muts_included <= 2+$mismatch);
				# $IS_flag = 0 if ($muts_rep >= $remain_spikein + 3 );
					
				# $NT_flag = 1 if ($remain_spikein == 0);
				# $NT_flag = 1 if ($count_ref >= 4);
				
				my $AL_flag=""; # Alternative alignments 
				if (int($array[1]) > 200){
					$AL_flag="AL";
				}
				
				if ($IS_flag == 1){
					if (exists $IS_reads{$array[0]}){
						$IS_reads{$array[0]}=$IS_reads{$array[0]}."\n".$AL_flag."IS-M".$info;
					}else{
						$IS_reads{$array[0]}=$AL_flag."IS-M".$info;
					}
				}elsif ($NT_flag == 1){
					# If dinulecotides detected, it immediately goes to native reads pool
					if (exists $NT_reads{$array[0]}){
						$NT_reads{$array[0]}=$NT_reads{$array[0]}."\n".$AL_flag."NT-M".$info;
					}else{
						$NT_reads{$array[0]}=$AL_flag."NT-M".$info;
					}
				}else{
					if (exists $SUS_reads{$array[0]}){
						$SUS_reads{$array[0]}=$SUS_reads{$array[0]}."\n".$AL_flag."SUS-M".$info;
					}else{
						$SUS_reads{$array[0]}=$AL_flag."SUS-M".$info;
					}
					$SUS_sam{$line."\n"}=1;
				}
			}else{
				# reads does not included spike-in mutation sites.
				if (exists $NT_reads{$array[0]}){
					$NT_reads{$array[0]}=$NT_reads{$array[0]}."\nNT-NS\t".$line;
				}else{
					$NT_reads{$array[0]}="NT-NS\t".$line;
				}
			}
		}else{
			# these reads are not mapped to provided IS reference.
			if (exists $NT_reads{$array[0]}){
				$NT_reads{$array[0]}=$NT_reads{$array[0]}."\nNT-UM\t".$line;
			}else{
				$NT_reads{$array[0]}="NT-UM\t".$line;
			}
		}
	}
}
close FH;

# remove doubt reads if (its paired reads) is already in IS/NT reads.
foreach (keys %IS_reads){
	if (exists $SUS_reads{$_}){
		$SUS_reads{$_} = $IS_reads{$_}."\n".$SUS_reads{$_};
		delete $IS_reads{$_};
	}
}

foreach (keys %NT_reads){
	if (exists $SUS_reads{$_}){
		$SUS_reads{$_} = $NT_reads{$_}."\n".$SUS_reads{$_};
		delete $NT_reads{$_};
	}
}

## Generate paired SAM FILE. ####
foreach (keys %SUS_sam){
    # Only doubted reads are go through the slow mode validation.
    print OFH_SAM_SUS $_;
}
close OFH_SAM_SUS;


# Double check suspected reads with human genome.
my $sus_filter_1 = $hq_header.".filter.1.fastq";
my $sus_filter_2 = $hq_header.".filter.2.fastq";
my %IS_reads_confirmed=();
my %SUS_reads_confirmed=();

#For suspected reads confirmation, first convert doubted sam to fastq inputs
system('java -jar '.$picard_file.' SamToFastq I='.$sam_file_sus." FASTQ=".$sus_filter_1." SECOND_END_FASTQ=".$sus_filter_2.
		" VALIDATION_STRINGENCY=LENIENT 2>>tmp/samtofastq.log");
if (-s $sus_filter_1 > 0){
	#mapping doubted reads to human reference
	open(FH,'bwa mem -h 1 -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_NT.' '.
		$sus_filter_1.' '.$sus_filter_2.' 2>>tmp/bwa_run_2.log|');
		
	while(<FH>){
		#ignore header information
		if ($_ !~ /^@/){
			my $flag_IS = 0;
			my $flag_NT = 0;
			
			my $line=$_;
			chomp($line);
			my @array=split("\t");
			my $md_NT = ();
			if (exists $SUS_reads{$array[0]} && (int($array[1]) < 200)){				
				#read MD string; only reads in doubt pool are examined.
				$md_NT = $md_NT."\t".$1 while ($line =~ /\tMD:Z:([\^\dACGT]+)\t/g); 
				
				my $ml_NT = 0;
				$ml_NT = $ml_NT + $1 while ($array[5] =~ /(\d+)M/g);
				
				my @array_IS = split("\t",$SUS_reads{$array[0]});
				my $ml_IS = 0;
				my $ml_IS = $ml_IS + $1 while ($array_IS[10] =~ /(\d+)M/g);
				
				# $flag_IS = 1 if ($md_NT =~ /[ATCG]0[ATCG]/ );
				$flag_NT = 1 if ($md_NT !~ /[ATCG]0[ATCG]/ && ($ml_NT > $map_score) && ($ml_NT > $ml_IS));
			}
			if ($flag_IS == 1){
				$IS_reads{$array[0]} = $SUS_reads{$array[0]};
				delete $SUS_reads{$array[0]};
			}
			if ($flag_NT == 1){
				$NT_reads{$array[0]} = $SUS_reads{$array[0]};
				my $NT_info = "HG_Mapping_Cigar: ".$array[5]." - HG_MD_String: ".$md_NT."\tSUS-";
				$NT_reads{$array[0]} =~ s/SUS-/$NT_info/g;
				delete $SUS_reads{$array[0]};
			}
		}
	}
	close FH;
}
## CLEAN UP 
system('rm -f '.$sus_filter_1);
system('rm -f '.$sus_filter_2);  

# Make output file
open(OFH_is,'>'.$hq_file_is);
open(OFH_sus,'>'.$hq_file_sus);
open(OFH_ntm,'>'.$hq_file_ntm);
open(OFH_ntns,'>'.$hq_file_ntns);
open(OFH_ntum,'>'.$hq_file_ntum);

# Native reads
my %NTM_reads = ();
my %NTNS_reads = ();
my %NTUM_reads = ();

foreach (keys %NT_reads){
	print OFH_ntm $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-M/ || $NT_reads{$_} =~ /SUS-/);
	print OFH_ntns $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-NS/ && $NT_reads{$_} !~ /NT-M/ && $NT_reads{$_} !~ /SUS-/);
	print OFH_ntum $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-UM/ && $NT_reads{$_} !~ /NT-M/ && $NT_reads{$_} !~ /NT-NS/ && $NT_reads{$_} !~ /SUS-/);
	$NTM_reads{$_} = $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-M/ || $NT_reads{$_} =~ /SUS-/);
	$NTNS_reads{$_} = $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-NS/ && $NT_reads{$_} !~ /NT-M/ && $NT_reads{$_} !~ /SUS-/);
	$NTUM_reads{$_} = $NT_reads{$_}."\n" if ($NT_reads{$_} =~ /NT-UM/ && $NT_reads{$_} !~ /NT-M/ && $NT_reads{$_} !~ /NT-NS/ && $NT_reads{$_} !~ /SUS-/);
	print $NT_reads{$_}."\n" if ($NT_reads{$_} !~ /NT-UM/ && $NT_reads{$_} !~ /NT-M/ && $NT_reads{$_} !~ /NT-NS/ && $NT_reads{$_} !~ /SUS-/);
}

close OFH_ntm;
close OFH_ntns;
close OFH_ntum;

# header/summary:
print OFH_is "Total Statistics:\n";
print OFH_is "1. Total Output Reads:\t".scalar(keys %reads_count_total)."\n";
print OFH_is "2. Total Mapped Reads:\t".scalar(keys %reads_count_bwa)."\n";
print OFH_is "3. IS-control Reads: ".scalar(keys %IS_reads)."\n";
print OFH_is "4. Suspected Reads: " .scalar(keys %SUS_reads)."\n";
print OFH_is "5. Mapped Native Reads: " .scalar(keys %NTM_reads)."\n";
print OFH_is "6. NonSite Native Reads: " .scalar(keys %NTNS_reads)."\n";
print OFH_is "7. UnMapped Native Reads: " .scalar(keys %NTUM_reads)."\n";
print OFH_is "8. Total Native Reads: " .scalar(keys %NT_reads)."\n";
# IS-control reads
foreach (keys %IS_reads){
	print OFH_is $IS_reads{$_}."\n";
}

# suspected reads
foreach (keys %SUS_reads){
	print OFH_sus $SUS_reads{$_}."\n";
}

close FH;   
close OFH_is;
close OFH_sus;

#output summary in the std-out
print  "Total Statistics:\n";
print  "1. Total Output Reads:\t".scalar(keys %reads_count_total)."\n";
print  "2. Total Mapped Reads:\t".scalar(keys %reads_count_bwa)."\n";
print  "3. IS-control Reads: ".scalar(keys %IS_reads)."\n";
print  "4. Suspected Reads: " .scalar(keys %SUS_reads)."\n";
print  "5. Mapped Native Reads: " .scalar(keys %NTM_reads)."\n";
print  "6. NonSite Native Reads: " .scalar(keys %NTNS_reads)."\n";
print  "7. UnMapped Native Reads: " .scalar(keys %NTUM_reads)."\n";
print  "8. Total Native Reads: " .scalar(keys %NT_reads)."\n";

print "Check the ".$hq_file_is." to see the details of captured spike-in reads.\n";

## generate fastq files
if ($type eq 'all' || $type eq 'gzip'){
    print "Done. Now generating output fastq/bam files ... \n";
    print "Warning: defined output files already exists; will overwrite! \n" if (-f $is_fastq_1);
    
    # initialize
    my $is_fastq_1_sam = $is_fastq_1;
    my $nt_fastq_1_sam = $nt_fastq_1;
    my $sus_fastq_1_sam = $sus_fastq_1;

    my $signal = 0;
    my $count_pair_1= 0;
    my $read_line_i = 0;
    
    # separating fastq file 1
    if($SAMPLE_1){  
        print ('Extracting Reads ... (Pair - 1)'."\n");
		if ($SAMPLE_1 =~ /\.gz$/){
			open(FH,'gzip -cd '.$SAMPLE_1.'|');
		}else{
			open(FH,'less '.$SAMPLE_1.'|');
		}    
    }elsif($BAM){
        print ('Extracting Reads from BAM ... '."\n");
        if ($keepBam){
            open(FH,$samtools_path.' view '.$BAM.'|');
            
            $is_fastq_1_sam =~ s/\.bam/\.sam/g;
            $nt_fastq_1_sam =~ s/\.bam/\.sam/g;
            $sus_fastq_1_sam =~ s/\.bam/\.sam/g;
			
        }else{
            open(FH,$samtools_path.' fastq '.$BAM.'|');
        }
    }
    if ($keepBam){
        open(OFH_spike,'>'.$is_fastq_1_sam);
        open(OFH_origin,'>'.$nt_fastq_1_sam);
        open(OFH_suspect,'>'.$sus_fastq_1_sam);

        open(FH_header, $samtools_path.' view -H '.$BAM.'|');
        while(<FH_header>){
            # output original sam header;
            print OFH_spike $_;
            print OFH_origin $_;
            print OFH_suspect $_;
        }
        close FH_header;
    }else{
        open(OFH_spike,'>'.$is_fastq_1);
        open(OFH_origin,'>'.$nt_fastq_1);
        open(OFH_suspect,'>'.$sus_fastq_1);
    }
    
	my $array = "";
    while(<FH>){
        #read data from the original fastq/bam files.
        my $line = $_;
        if($keepBam){
            my @array = split("\t",$line);
            print OFH_spike $line if (exists $IS_reads{$array[0]});
            print OFH_origin $line if (exists $NT_reads{$array[0]});
            print OFH_suspect $line if (exists $SUS_reads{$array[0]});
        }else{
            if ($read_line_i % 4 ==0){ # only test on the head line
                # Casava 1.8 format || Illumina reads ends with /1 or /2
                if (/^@(\S+)\s(\S+)/ || /^@(\S+)\/[12]$/){
                    $array=$1;
                }elsif(/^@(.+)/){
                    $array=$1;
                    chomp($array);
                }
            }
            $read_line_i = $read_line_i + 1 ;
            print OFH_spike $line if (exists $IS_reads{$array});
            print OFH_origin $line if (exists $NT_reads{$array});
            print OFH_suspect $line if (exists $SUS_reads{$array});
        }
    }
    close FH;
    close OFH_spike;
    close OFH_origin;
    close OFH_suspect;
        
    # If paired-end reads;
    if($SAMPLE_2 ne 'empty'){
        # separating fastq file 2
        print ('Extracting Reads... (Pair - 2)'."\n");
        if ($SAMPLE_2 =~ /\.gz$/){
			open(FH,'gzip -cd '.$SAMPLE_2.'|');
		}else{
			open(FH,'less '.$SAMPLE_2.'|');
		} 
        open(OFH_spike,'>'.$is_fastq_2)  ;
        open(OFH_origin,'>'.$nt_fastq_2);
        open(OFH_suspect,'>'.$sus_fastq_2);
        
        $signal = 0;
        $read_line_i = 0;
		my $array = "";
        while(<FH>){
            my $line = $_;
            if ($read_line_i % 4 ==0){ # head line
                # Casava 1.8 format || Illumina reads ends with /1 or /2
                if (/^@(\S+)\s(\S+)/ || /^@(\S+)\/[12]$/){
                    $array=$1;
                }elsif(/^@(.+)/){
                    $array=$1;
                    chomp($array);
                }
            }
            $read_line_i = $read_line_i + 1 ;
            print OFH_spike $line if (exists $IS_reads{$array});
            print OFH_origin $line if (exists $NT_reads{$array});
            print OFH_suspect $line if (exists $SUS_reads{$array});
        }
        close FH;
        close OFH_spike;
        close OFH_origin;
        close OFH_suspect;
    }
    
    if($keepBam){
        # Convert sam to bam file.
        system($samtools_path.' view -b '.$is_fastq_1_sam.' >'.$is_fastq_1);
        system($samtools_path.' view -b '.$nt_fastq_1_sam.' >'.$nt_fastq_1);
        system($samtools_path.' view -b '.$sus_fastq_1_sam.' >'.$sus_fastq_1);
        # Remove original sam file.
        system('rm -f '.$is_fastq_1_sam);
        system('rm -f '.$nt_fastq_1_sam);
        system('rm -f '.$sus_fastq_1_sam);
    }
}   
# End of the processing.}
print ('Done! In total '.scalar(keys %IS_reads) ." Spike-in reads have been separated.\n");

######   Generating Gzip files ######
if ($type eq 'gzip'){
    print ('Gzip output files.. (may take long time...)'."\n");
    system('gzip -1 '.$is_fastq_1);
    system('gzip -1 '.$nt_fastq_1);
    system('gzip -1 '.$sus_fastq_1);
    if($SAMPLE_2){
        system('gzip -1 '.$is_fastq_2);
        system('gzip -1 '.$nt_fastq_2);
        system('gzip -1 '.$sus_fastq_2);
    }
    print ('Program finished normally.'."\n");
}

## Main Program End

## Subroutines for Test mode                 
sub get_test_data {
    my ($test_mode) = @_;
    my $SAMPLE_1 = '';
    my $SAMPLE_2 = 'empty';
    my $BAM = '';
    my $hq_header = '';
    print("In Test Mode\n");
    
    if ($test_mode eq 'test_pos') {
        # test Spike-in Sample
        $SAMPLE_1='data/test_pos_2_1.fastq.gz';
        $SAMPLE_2='data/test_pos_2_2.fastq.gz';
        
        $hq_header = $SAMPLE_1; 
        $hq_header =~ s/_1\.fastq\.gz//g;
        $hq_header =~ s/.*\///g;
    }
    if ($test_mode eq 'test_neg') {     
        # test 'Negative' Sample
        $SAMPLE_1='data/test_neg_2_1.fastq.gz'; 
        $SAMPLE_2='data/test_neg_2_2.fastq.gz';

        $hq_header = $SAMPLE_1; 
        $hq_header =~ s/_1\.fastq\.gz//g;
        $hq_header =~ s/.*\///g;
    }
    if ($test_mode eq 'test_bam') {     
        # test 'Positive' Sample in BAM format
        $BAM='data/test_pos.bam'; 

        $hq_header = $BAM;
        $hq_header =~ s/\.bam/_bam/g;   
        $hq_header =~ s/.*\///g;
    }
    return ($SAMPLE_1, $SAMPLE_2, $BAM, $hq_header);
}

sub prepare_ref {
    my ($ref_prefix) = @_;
	my $IS_prefix = $ref_prefix.'.IS';
	my $NT_prefix = $ref_prefix.'.NT';
	print ("customized spike-in reference are used. (".$ref_prefix.") Check required files ... \n");
    # Create the $mut_ref file of new spike-in reference if not existed.
    my $mut_REF=$IS_prefix.'.vcf';
    my $mut_bed=$IS_prefix.'.bed';
    if (! -f $mut_REF){
        print ("No Mutation position file detected, try generating them automatically... \n");
        open(FH,$IS_prefix.'.fasta');
        open(OFH,'>', $mut_REF);
        open(OFH_bed,'>', $mut_bed);
        # Initialize
        my $curr_gene_name = '';
        my $i = 0;
        
        while( <FH> ){
			my $line = $_; 
            if (/^>(\S+)/){
                # reset position to 0 when a new gene/chrom started.
                $curr_gene_name = $1;
                $i = 0;
			}else {
                my $string = $_;
                chomp($string);
                foreach (split //, $string){
                    $i = $i + 1;
                    if ($_ =~ /[ACGT]/){
                        print OFH $curr_gene_name."\t".$i."\t".$_."\n";
                        print OFH_bed $curr_gene_name.":".$i."-".$i."\n";
                    }
                }
			}
        }
        
        close FH;
        close OFH;
        close OFH_bed;
        
        if (-f $mut_REF){
            print("Success.\n");
        }else{
            print("Failed... please contact the developer.\n");
            exit 0;
        }
    }else{
        print ("Mutation position file found!\n");
    }
	
	# Create BWA INDEX if not existed.
    if (! -f $IS_prefix.'.fasta.bwt'){
        print ("BWA index files are not detected, try generating them automatically... \n");
        system('bwa index '.$IS_prefix.'.fasta '); 
        if (-f $IS_prefix.'.fasta.bwt'){
            print ("Successfully generate BWA index. \n");
        }else{
            print ("Failed, please check BWA availability ... \n");
            exit 0;
        }
    }else{
        print ("BWA INDEX found!\n");
    }
	
	if (! -f $NT_prefix.'.fasta.bwt'){
        print ("BWA index files are not detected, try generating them automatically... \n");
        system('bwa index '.$NT_prefix.'.fasta '); 
        if (-f $NT_prefix.'.fasta.bwt'){
            print ("Successfully generate BWA index. \n");
        }else{
            print ("Failed, please check BWA availability ... \n");
            exit 0;
        }
    }else{
        print ("BWA INDEX found!\n");
    }
    
    return $mut_REF;
}