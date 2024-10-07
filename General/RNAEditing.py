#!/usr/bin/python
import gzip
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-r", "--rna", help = "RNA VCF file", required=True)
parser.add_argument("-d", "--dna", help = "DNA VCF file", required=True)
parser.add_argument("-o", "--output", help = "output file for the results (required)", required=True)
args = parser.parse_args()

rna=gzip.open(args.rna, "r")
dna=gzip.open(args.dna, "r")
out=open(args.output, "w")
D={}

dp_cutoff=20
dp_id=4

print "building dictionary from RNA VCF"
for line in rna:
        if line.startswith("#"):
                continue
        else:
                tokens=line.split("\t")
                info=tokens[7].split(";")
                if info[3].split("=")[0]=="DP":
                        dp_id=3
                elif info[4].split("=")[0]=="DP":
                        dp_id=4
                dp=float(info[dp_id].split("=")[1])
                alt=tokens[4]
                ref=tokens[3]
                if dp >= dp_cutoff and ((ref=="A" and alt=="G") or (ref=="T" and alt=="C")):
                        alt_dp=float(tokens[9].split(":")[1].split(",")[1])
                        edit_ratio=alt_dp/dp
                        key=tokens[0].strip("chr").strip("Chr")+"_"+tokens[1]
                        #print line
                        #print "alt_dp is"+ str(edit_ratio)
                        value=alt+"_"+str(edit_ratio)
	                D[key]=value
                
print "comparing to dictionary with second VCF"
for line in dna:
        if line.startswith("#"):
                continue
        else:
                tokens=line.split("\t")
                key=tokens[0].strip("chr").strip("Chr")+"_"+tokens[1]
                gts=tokens[9].split(":")
                gt=gts[0]
                alt=tokens[4].split(",")[0]
                dp=gts[1]
                ref=tokens[3]
                alt=tokens[4]
                if (ref=="A"  or ref =="T") and alt.startswith("<") and gt=="0/0" and key in D:
                        if int(dp) >= dp_cutoff:
                                a1=tokens[3]
                                value=D[key].split("_")
                                a2=value[0]
                                edit_ratio=value[1]
                                if (a1=="A" and a2=="G") or (a1=="T" and a2=="C"):
                                        print "EDITING: "+tokens[0]+"\t"+tokens[1]+"\t"+tokens[1]+"\t" + a1+"\t"+a2 +"\t"+dp+"\t"+edit_ratio
                                        out.write(tokens[0]+"\t"+tokens[1]+"\t"+tokens[1]+"\t" + a1+"\t"+a2+"\t"+dp+"\t"+edit_ratio+"\n")

rna.close()
dna.close()
out.close

