## Trans-Bias

Transcription Bias

## Installation

Trans-Bias was implemented using a combination of R and Python3. 

```
sudo apt update
sudo apt install pandoc
```
If you don't have pip installed:

``` 
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py 
sudo python get-pip.py
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

Vcf, bam, gtf and Tdp arguments are mandatory. Users must prepare one tab delimited vcf file with sample headers the same as bam filenames. For example, a vcf file with sample column name NAP1, NAP2, NAP3 must have corresponding bam files named NAP1.bam, NAP2.bam, NAP3.bam. Tdp is the Total depth value given by the user to subset the vcf. Gtf file must correspond to the genome reference used for alignment. 

You can download example bam and gtf files:

```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1.gtf
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/Bam
```
