import argparse
import glob
import os
import shutil
import json
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import time
start = time.time()

def int_pos(value):
    i = int(value)
    if i <= 0:
        raise argparse.ArgumentTypeError("%s is not a positive integer value" % value)
    return i
	
# arguments parsing
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
parser.add_argument('Tdp', help = "Input the total depth threshold", type=int_pos)
parser.add_argument("-v", "--vcf", help = "Path to the vcf file (required)", required=True)
parser.add_argument("-b", "--bam", help='Path to the bam file (required)', required = True)
parser.add_argument("-g", "--gtf", help= 'Path to the gtf file (required).', required = True)

args = parser.parse_args()
    
num_sample=0
bed=[]

currentDir = os.getcwd()
bedDir = currentDir + '/bedfiles/'
if not os.path.exists(bedDir):
	os.mkdir(bedDir)
    
pileupDir= currentDir + '/pileups/'
if not os.path.exists(pileupDir):
	os.mkdir(pileupDir)
    
outputDir= currentDir + '/results/'
if not os.path.exists(outputDir):
	os.mkdir(outputDir)


#loading vcf files then output bed
def step1(vcf, bedD, vcfdepth):
    print("Starting Step 1")
    print("\tVcf File Location: " + vcf)
    print("\tBedfiles Directory: " + bedD)

    try:
        input = open(vcf, "r")
    except (OSError, IOError) as e:
        print("####ERROR:can not load vcf file " + vcf)
        raise Exception(e)

    for line in input:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            samples=line.split("\t")
            num_sample=len(samples)-9
            print("\t\tThere are "+ str(num_sample)+ " samples in " + vcf + "\n")
            for i in range(9, len(samples)):
                filename=samples[i]+".bed"
                output_bed=open(bedD+filename, "w")
                bed.append(output_bed)
            continue
      
        words=line.strip("\n").split("\t")
        for i in range(9, len(words)):
            tokens=words[i].split(":")
            dp=0
            if len(tokens) > 1:
                dp_field=tokens[1].split(",")
                if len(dp_field) >1:
                    dp=int(dp_field[0])+int(dp_field[1])
            if tokens[0] =="0/0" or tokens[0]=="./." or dp < vcfdepth:
                continue
            bed[i-9].write(words[0]+"\t"+words[1]+"\t"+words[3]+"\n")

#############################
if os.path.exists(args.vcf):
    step1(args.vcf, bedDir, args.Tdp)
else:
    print("Vcf file not found")
    exit()  
