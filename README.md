#Transcription Bias

## Installation

Trans-Bias was implemented using a combination of R and Python3. 

Install pandoc as root
```
sudo apt update
sudo apt install pandoc
```

Ensure Python version >=3.8 then install:
```
pip install biopython
```

Ensure R version is >=3.6 then install:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")

install.packages(c("rmarkdown", "knitr", "jsonlite", "ggplot2", "reshape2"))
```

## Usage:

Download from Github
```
git clone https://github.com/MaryOlufunmilola/Trans-Bias

usage: python3 TBias.py [-h] -v [VCF] -b [BAM] -g [GTF] Tdp

positional arguments:
  Tdp                Input the total depth threshold

optional arguments:
  -h,           --help              show this help message and exit
  -v [VCF],     --vcf [VCF]         vcf file (required)
  -b [BAMDIR],  --bamDir [BAMDIR]   Full absolute Path to the bam file(s) (required)
  -g [GTF],     --gtf [GTF]         gtf file (required).
```

For example:
```
cd Trans-Bias
python3 TBias.py 20 -v example.vcf -b ~/Trans-Bias/Bam/ -g chr1.gtf
```

Vcf, bam, gtf and Tdp arguments are mandatory. Users must prepare one tab delimited vcf file with sample headers the same as bam filenames. For example, a vcf file with sample column name NAP1, NAP2, NAP3 must have corresponding bam files named NAP1.bam, NAP2.bam, NAP3.bam. Tdp is to limit genomic sites to those with read Total depth > value given by the user. Gtf file must correspond to the genome reference used for alignment. 

You can download example bam and gtf files:
```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1.gtf
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/Bam
```

## How it works:
STEP1: Generate bed file for each sample in vcf file (remove genotype unphased call '0/0', missing allele './.' and Total Depth  < user defined at each position from vcf samples)

STEP2: Generate pileup from Bam and bed files using Rsamtools

STEP3: Converts pileup output into a table of nucleotide frequencies at each SNP position based on strand orientation (A,T,C,G,a,c,g,t)

STEP4: Gene and strand annotatation using reference gtf file (**Genename is not used for further steps at the moment?**)

STEP5: Get alternate allele (usually nucleotides with maximum count value or second maximum if maximum is same as ref) for transcription and non transcription

STEP6: Sum values in rows with same ref to alt allelle for transcription and non transcription

STEP7: Generate figures