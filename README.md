## Trans-Bias

Transcription Bias analysis

### Installation

Trans-Bias was implemented using a combination of R and Python3 codes. 
If you don't have pip installed:

``` 
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py 
sudo python get-pip.py
```

Check if Python3.7 is installed in Ubuntu.

```
python -version
```

If Python version is lower than 3.8 or not installed, run the commands below sequentially:

```
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.8
```

Install R packages in >=R3.6

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")

install.packages(c("rmarkdown", "knitr", "jsonlite", "ggplot2", "reshape2"))
```

### Usage:

Download from Github

```
git clone https://github.com/MaryOlufunmilola/Trans-Bias
cd TBias 

usage: TBias_13456.py [-h] -v [VCF] -b [BAM] -g [GTF] Tdp

positional arguments:
  Tdp                Input the total depth threshold

optional arguments:
  -h, --help         show this help message and exit
  -v [VCF], --vcf [VCF]  vcf file (required)
  -b [BAMDIR], --bamDir [BAMDIR]  Full Path to the bam file(s) (required)
  -g [GTF], --gtf [GTF]  gtf file (required).
```

For example:
```
nohup python3 TBias.py 20 \
    -v example.vcf \
    -b Trans-Bias-Bam/ \
    -g chr1.gtf &
```

You can download example bam, vcf and gtf files:

```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/chr1.gtf
```