#!/usr/bin/python
import sys

ad={}
df1=set()
df2=set()
header1 = "CHR"+"\t"+"START"+"\t"+"END"+"\t"+"GENENAME"+"\t"+"GENEBIOTYPE"
header2 = "CHR"+"\t"+"START"+"\t"+"END"+"\t"+"GENENAME"
prefixes = ("gene_name", "gene_biotype")
not_needed = ['misc_RNA', 'scaRNA', 'TEC', 'Mt_tRNA', 'sense_overlapping', 'sRNA', 'sense_intronic', 'non_coding', 'ribozyme', 'processed_transcript', 'Mt_rRNA']

with open(sys.argv[1], 'r') as f1:
    for line in f1:
	if line.startswith("#!"):
            continue
	fs = line.split("\t")
        chrom=fs[0]
        region=fs[2]
        st=fs[3]
        en=fs[4]
        wt = fs[8].split("; ")
        genes = [x for x in wt if x.startswith(prefixes)]
        gs=genes[0].split(" ")
        bt=genes[1].split(" ")
        gene=gs[1].replace('"', '')
        biotype=bt[1].replace('"', '').replace(';\n', '')
        ad[chrom+'_'+st+'_'+en] = [region,gene,biotype]
        if region not in df1:
            df1.add(region)
            d1=list(df1)
            d1 = [e for e in d1 if e not in ('gene', 'transcript', 'exon')]
            
        if biotype not in df2:
            df2.add(biotype)
            d2=list(df2)
            d2= list(set(["pseudogene" if x.endswith("pseudogene") else x for x in d2]))
            d2[:]= (y for y in d2 if not y.endswith("_gene"))
            d2 = [e for e in d2 if e not in (not_needed)]
                        

for i in range(0,len(d1)):
    filenames1=d1[i]+".bed"
    with open(filenames1, 'w') as f3:
        f3.write(header1+"\n")
        for k,v in ad.items():
            if str(d1[i]) == str(v[0]):
                tw=k.split("_")
                f3.write(str(tw[0])+"\t"+str(tw[1])+"\t"+str(tw[2])+"\t"+str(v[1])+"\t"+str(v[2])+"\n")

for i in range(0,len(d2)):
    filenames2=d2[i]+"_exon.bed"
    with open(filenames2, 'w') as f4:
        f4.write(header2+"\n")
        for k,v in ad.items():
            if str(d2[i]) == str(v[2]):
                tw=k.split("_")
                f4.write(str(tw[0])+"\t"+str(tw[1])+"\t"+str(tw[2])+"\t"+str(v[1])+"\n")
