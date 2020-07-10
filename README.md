## Trans-Bias

Transcription analysis of Aligned reads. Trans-Bias depends on R and Python 3.

## Installation

If you don't have pip installed:

``` 
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py 
sudo python get-pip.py
```

Check if Python3.7 is installed in Ubuntu.

```
python -version
```

If Python version is lower than 3.7 or not installed, run the commands below sequentially:

```
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.7
```

Install R packages in R

```
install.packages(c("optparse", "rmarkdown", "knitr", "pscl", "gridExtra"))
```

## Usage

Download from Github

```
git clone https://github.com/MaryOlufunmilola/Trans-Bias
cd TBias 

usage: TBias_13456.py [-h] -v [VCF] -b [BAM] -g [GTF] Tdp

positional arguments:
  Tdp                Input the total depth threshold

optional arguments:
  -h, --help         show this help message and exit
  -v [VCF], --vcf [VCF]  Path to the vcf file (required)
  -b [BAM], --bam [BAM]  Path to the bam file (required)
  -g [GTF], --gtf [GTF]  Path to the gtf file (required).
```

For example:
```
nohup python3 TBias_13456.py 20 \
    -v chr1.vcf \
    -b chr1_Bam/ \
    -g chr1.gtf
```

You can download example bam, vcf and gtf files:

```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/example1.bam
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/example1.bam.bai
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/example.vcf
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/chr1.gtf
```
