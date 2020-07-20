## Trans-Bias

Transcription analysis of Aligned reads. Trans-Bias depends on R and Python 3.

## Installation

If Python version is lower than 3.8 or not installed, run the commands below sequentially:

```
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.8
```

If you don't have pip installed:

``` 
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py 
sudo python3.8 get-pip.py
python3.8 -m pip install biopython
```

Install R packages in >=R3.4

```
BiocManager::install("Rsamtools")
install.packages(c("Rsamtools", "reshape2","ggplot2","jsonlite","rmarkdown"))
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
  -v [VCF], --vcf [VCF]  vcf file (required)
  -b [BAMDIR], --bamDir [BAMDIR]  Full Path to the bam file(s) (required)
  -g [GTF], --gtf [GTF]  gtf file (required).
```

For example:
```
nohup python3 TBias_13456.py 20 \
    -v example.vcf \
    -b chr1_Bam/ \
    -g chr1.gtf
```

You can download example gtf file:

```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans-Bias/chr1.gtf
```
