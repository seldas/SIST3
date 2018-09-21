# SIST: Separate Accugenomics Spike-In reads
## Version 3.0.1
09/21/2018: version updated to 3.0.1  
 -- added one option named -m (-maxmismatch) which is the threshold of how many IS dinulecotide position are allowed to mismatched in one sequence, default is 1.  
 -- Changed the default reference files. In version 3 the human reference genome is not needed as we didn't see much improvement of performance to use the whole genome. Meanwhile, the native sequence of dinucleotides are used as the native reference genome.  
 -- Picard are required by the tool, as well as samtools and bwa. the location of picard could be customized as option (-picard).


## Introduction

Spike-In Separation Toolbox (SIST) is developed for separate Accugenomics sequencing spike-in controls (reads） from original sample reads. The separation is majorly based on the unique characteristics of spike-in controls, which contains several di-nucleotides inside the target region.  

This toolbox has been specificially designed and used in *Sequencing Quality Control 2* project - Working Group II for the purpose of target sequencing quality control. To read more details about SEQC2 project, please visit [here](https://www.fda.gov/ScienceResearch/BioinformaticsTools/MicroarrayQualityControlProject/ucm507935.htm)   

For version 3, we allowed to output the third bin, for suspected reads. Suspected reads are those reads neither mapped to IS reference nor to NT reference, which basically included more errors in the dinucleotides position.

## Tool Manual (Help) Page
use `perl SIST3.pl -help` to get help information.

#### **Note**
During to privacy concern, currently the reference of Accugenomics Spike-in Control are not publicly released. This tool **CANNOT** be ran successfully without the specific refernece genome. If you are project related person, please contact [Leihong](mailto:leihong.wu@fda.hhs.gov) for further assistance.  

##Usage:
      SIST3.pl [-fastq_1=<fastq.gz> -fastq_2=<fastq.gz>] [-bam=<bam>] -O=<output_prefix> [options]
    
    Required inputs:
      -fastq_1  paired-end reads 1 (required if no bam input)
      -fastq_2  paired-end reads 2 (required if paired-end reads)
    
      -bam  input file as bam format (required if no fastq input)
    
      -O    Output file header 
    
    Common options:
      -ref [prefix]: prefix of reference, (defalt: "Refs/ISREF");
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
      -picard : picard path (default: Refs/picard.jar);
      
    -help:  Prints out this helpful message

## Example 
 
`perl SIST3.pl -f1=data/test_1.fastq.gz -f2=data/test_2.fastq.gz -O=output/test`


## Key algorithms
We used FLAGs to determine whether a read is IS or NT reas based on following criteria:  

(1) $muts_included: IS position (dinucleotides) in the test reads.  
(2) $count_ref:  How many REF (NT) occurred in these positions.  
(3) $remain_spikein: How many ALT (IS) occurred in these positions.
(4) $mismatch: customized threshold of option -m;

$IS_flag = 1 if ($remain_spikein >= $muts_included-$mismatch && $muts_included > 2+$mismatch);  
$IS_flag = 1 if ($remain_spikein == $muts_included && $muts_included >=2 && $muts_included <= 2+$mismatch);  

$NT_flag = 1 if ($count_ref >= $muts_included-$mismatch  && $muts_included > 2+$mismatch );  
$NT_flag = 1 if ($count_ref == $muts_included  && $muts_included >=2 && $muts_included <= 2+$mismatch);  

if the read had both IS_flag and NT_flag ==0, it will go to the Native reads. Or, if the read had $muts_included==0 or didn't mapped to the IS sequence, it will also be considered as Native reads.  

Then, we used the NT reference for further fine-tuning. as if the suspected reads mapped much better to NT reference than IS reference and didn't have any dinucleotides, it will also be considered as NT reads.
  

				

