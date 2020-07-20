import argparse
import glob
import os
import shutil
import json
from argparse import RawTextHelpFormatter
from Bio import SeqIO

def int_pos(value):
    i = int(value)
    if i <= 0:
        raise argparse.ArgumentTypeError("%s is not a positive integer value" % value)
    return i
	
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter) 
parser.add_argument('Tdp', help = "Input the total depth threshold", type=int_pos)
parser.add_argument("-v", "--vcf", help = "vcf file (required)", required=True)
parser.add_argument("-b", "--bamDir", help='Full Path to the bam file(s) (required)', required = True)
parser.add_argument("-g", "--gtf", help= 'gtf file (required).', required = True)

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

def step1(vcf, bedD, vcfdepth):
    print("Starting Step 1")
    print("\tVcf File Location: " + vcf)
    print("\tBedfiles Directory: " + bedD)

    try:
        input = open(vcf, "r")
    except (OSError, IOError) as e:
        print("ERROR:can not load vcf file " + vcf)
        raise Exception(e)

    for line in input:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            samples=line.rstrip('\n').split("\t")
            num_sample=len(samples)-9
            print("\tThere are "+ str(num_sample)+ " samples in " + vcf + "\n")
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

def step2(step, bam, bed, pileup):
	print("Starting Step 2")
	print("\tBamfile Directory: " + bam)
	print("\tBedfile Directory: " + bed)
	print("\tOutput Pileup Directory: " + pileup)
	os.system('/root/R_alternatives/3.6.2/bin/Rscript TBias_278.r ' + step + ' ' + bam + ' ' + bed + ' '+ pileup)

def step3(pileup):
	print("Starting Step 3")
	print("\tPileup Directory: " + pileup)

	gg=glob.glob(pileup+"*.pileup")
	dict2 ={}
	for x in gg:
		plist2 = []
		f=open(x, "r")
		next(f)
		l = f.readline()
		words=l.strip("\n").split("\t")
		key0=words[0]+"_"+words[1]+"_"+words[6]

		A,C,T,G,a,c,t,g = [0 for _ in range(8)]

		if words[3]=="A" and words[2]=="+":
			A=words[4]
		elif words[3]=="A" and words[2]=="-":
			a=words[4]
		elif words[3]=="T" and words[2]=="+":
			T=words[4]
		elif words[3]=="T" and words[2]=="-":
			t=words[4]
		elif words[3]=="C" and words[2]=="+":
			C=words[4]
		elif words[3]=="C" and words[2]=="-":
			c=words[4]
		elif words[3]=="G" and words[2]=="+":
			G=words[4]
		elif words[3]=="G" and words[2]=="-":
			g=words[4]

		for line in f:
			words=line.strip("\n").split("\t")
			key1=words[0]+"_"+words[1]+"_"+words[6]
    
			if key1==key0:
				if words[3]=="A" and words[2]=="+":
					A=words[4]
				elif words[3]=="A" and words[2]=="-":
					a=words[4]
				elif words[3]=="T" and words[2]=="+":
					T=words[4]
				elif words[3]=="T" and words[2]=="-":
					t=words[4]
				elif words[3]=="C" and words[2]=="+":
					C=words[4]
				elif words[3]=="C" and words[2]=="-":
					c=words[4]
				elif words[3]=="G" and words[2]=="+":
					G=words[4]
				elif words[3]=="G" and words[2]=="-":
					g=words[4]
                   
			else:
				tokens=key0.split("_")
				plist2.append(tokens[0]+"\t"+ tokens[1]+"\t"+ tokens[2]+"\t"+str(A)+"\t"+str(T)+"\t"+str(C)+"\t"+str(G)+"\t"+str(a)+"\t"+str(t)+"\t"+str(c)+"\t"+str(g)+'\n')
				A,C,T,G,a,c,t,g = [0 for _ in range(8)]

				key0=key1
				if words[3]=="A" and words[2]=="+":
					A=words[4]
				elif words[3]=="A" and words[2]=="-":
					a=words[4]
				elif words[3]=="T" and words[2]=="+":
					T=words[4]
				elif words[3]=="T" and words[2]=="-":
					t=words[4]
				elif words[3]=="C" and words[2]=="+":
					C=words[4]
				elif words[3]=="C" and words[2]=="-":
					c=words[4]
				elif words[3]=="G" and words[2]=="+":
					G=words[4]
				elif words[3]=="G" and words[2]=="-":
					g=words[4]

		plist2.append(words[0]+"\t"+ words[1]+"\t"+ words[6]+"\t"+str(A)+"\t"+str(T)+"\t"+str(C)+"\t"+str(G)+"\t"+str(a)+"\t"+str(t)+"\t"+str(c)+"\t"+str(g) + '\n')
		dict2[x] = plist2
	return dict2

def step4(gtf):
    print("Starting Step 4")
    print("\tGtf file location: " + gtf)

    dict4 = dict()
    for x in dict2:
        plist3 = []
        dict3 = dict()
        
        for line in dict2[x]:
            b = line.split('\t')
            dict3[b[0]+'_'+b[1]]=line.strip("\n")
 
        with open(gtf, 'r') as f:
            for line in f:
                line= line.rstrip(os.linesep)
                if not line.startswith('#!'):
                    lss=line.split('\t')
                    if lss[2] == "gene":
                        for i in range(int(lss[3]), int(lss[4])+1):
                            tokens=lss[8].split(";")
                            gene_name=tokens[1].split('\"')
                            gene_name=gene_name[1]
                            key=lss[0].lstrip(' ')+"_"+str(i)
                            if key in dict3.keys():
                                plist3.append(dict3[key]+"\t"+gene_name+"\t"+lss[6]+"\n")
                                
        dict4[x] = plist3
    return dict4

def step5():
    print("Starting Step 5")

    dict5 = dict()
    for x in dict4:
        plist4 = []
        tr,ntr,Aa,Cc,Tt,Gg = [0 for _ in range(6)]
        ref_alt=""
        list1 = []
        prev_line=''
        for line in dict4[x]:
            if line != prev_line:
                w=line.strip("\n").split("\t")
                key1=w[11]+"_"+w[12]
                
                Aa = w[3]+w[7]
                Tt = w[4]+w[8]
                Cc = w[5]+w[9]
                Gg = w[6]+w[10]
                list1 = [Aa,Tt,Cc,Gg]
                
                maxNum = max(list1) 
                maxIndex = list1.index(max(list1))
                secondNum=min(list1) 
                secondIndex = list1.index(min(list1))
                
                for index, number in enumerate(list1):
                    if index == maxIndex: continue
                    if number > secondNum:
                        secondNum = number
                        secondIndex = index
                
                if w[12] == "+":
                    if secondIndex == 0:
                        alt1 = "A"
                        tr1 = w[3]
                        ntr1 = w[7]
                    elif secondIndex == 1:
                        alt1 = "T"
                        tr1 = w[4]
                        ntr1 = w[8]
                    elif secondIndex == 2:
                        alt1 = "C"
                        tr1 = w[5]
                        ntr1 = w[9]
                    elif secondIndex == 3:
                        alt1 = "G"
                        tr1 = w[6]
                        ntr1 = w[10]
                    
                    if maxIndex == 0:
                        alt = "A"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[3]
                            ntr=w[7]
                    
                    elif maxIndex == 1:
                        alt = "T"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[4]
                            ntr=w[8]
                    
                    elif maxIndex == 2:
                        alt = "C"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[5]
                            ntr=w[9]
                    
                    elif maxIndex == 3:
                        alt = "G"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[3]
                            ntr=w[7]
                    

                else:
                    if secondIndex == 0:
                        alt1 = "A"
                        tr1 = w[7]
                        ntr1 = w[3]
                    elif secondIndex == 1:
                        alt1 = "T"
                        tr1 = w[8]
                        ntr1 = w[4]
                    elif secondIndex == 2:
                        alt1 = "C"
                        tr1 = w[9]
                        ntr1 = w[5]
                    elif secondIndex == 3:
                        alt1 = "G"
                        tr1 = w[10]
                        ntr1 = w[6]
                    
                    if maxIndex == 0:
                        alt = "A"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[7]
                            ntr=w[3]
                    
                    elif maxIndex == 1:
                        alt = "T"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[8]
                            ntr=w[4]
                    
                    elif maxIndex == 2:
                        alt = "C"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[9]
                            ntr=w[5]
             
                    elif maxIndex == 3:
                        alt = "G"
                        if w[2] == alt:
                            ref_alt=w[2]+"->"+alt1
                            tr=tr1
                            ntr=ntr1
                        else:
                            ref_alt=w[2]+"->"+alt
                            tr=w[10]
                            ntr=w[6]
                prev_line = line
                plist4.append(str(w[0])+"\t"+str(w[1])+"\t"+str(ref_alt)+"\t"+str(tr)+"\t"+str(ntr)+"\t"+str(w[11])+"\t"+ str(w[12])+'\n')
        
        dict5[x] = plist4
    return dict5
    
def step6():
	print("Starting Step 6")
	dict6 = {}
	for x in dict5:
		dict6[x] = {}
		for line in dict5[x]:
			v=line.strip("\n").split("\t")
			values = [int(i) for i in v[3:5]]
			if v[2] not in dict6[x]:
				dict6[x][v[2]] = values
			else:
				dict6[x][v[2]] = [sum(k) for k in zip(dict6[x][v[2]],values)] 
		print(dict6[x])
	return dict6

def step7(step, output):
	print("Starting Step 7")
	print("Creating markdown summary")
	print("\tOutput Directory: " + output)
	os.system('/usr/bin/Rscript TBias_278.r ' + step + ' ' + output)

#############################
if os.path.exists(args.vcf):
    step1(args.vcf, bedDir, args.Tdp)
else:
    print("Vcf file not found")
    exit()  

if not os.listdir(args.bamDir) :
    print("Directory is empty")
    exit()
else:    
    step2("2", args.bamDir, bedDir, pileupDir)
	
if any(File.endswith(".pileup") for File in os.listdir(pileupDir)):
    dict2 = step3(pileupDir)
    shutil.rmtree(bedDir)
else:
    print("Pileup files are not in this directory")
    exit()
    
if dict2:
    dict4 = step4(args.gtf)
    shutil.rmtree(pileupDir)
	
else:
    print("Gtf file not in this directory, Cannot move on to step 4")
    exit()
    
if dict4:
    dict5 = step5()
    dict2.clear()
else:
    print("Cannot move on to step 5, Error in previous step")
    exit()
    
if dict5:
    dict6 = step6()
    with open('dict5.json', 'w') as f1:
        json.dump(dict5, f1)
    with open('dict6.json', 'w') as f2:
        json.dump(dict6, f2)
else:
    print("Cannot move to Step 6, Error in previous step")
    exit()

if dict6:
    step7("7", outputDir)
    dict5.clear()
    dict6.clear()
else:
    print("Cannot move to Step 7, Error in previous step")
    exit()
    