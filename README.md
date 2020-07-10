## Trans-Bias

Transcription analysis of Aligned reads. Trans-Bias depends on R and Python 3.

## Installation

Download from Github

```
git clone https://github.com/MaryOlufunmilola/Trans-Bias
```

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

nohup python3 TBias_steps.py depth_value -v vcf_file -b bam_file -g gtf_file &

For example:
```
nohup python3 TBias_steps.py 20 \
    -v chr1.vcf \
    -b chr1.bam \
    -g chr1.gtf
```

You can download our example bam, vcf and gtf files:

```
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1_Bam/No_Treatment_5.bam
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1_Bam/No_Treatment_5.bam.bai
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1.vcf
wget http://www.innovebioinfo.com/Sequencing_Analysis/Trans_Bias/chr1.gtf
```
