#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
step=strtoi(args[1])

if (step == 2) {
    Bampath= args[2]
    Bedpath= args[3]
    Outpath= args[4]
    setwd(Outpath)
	
	library(Rsamtools)

    files_bam = list.files(path = Bampath, pattern="\\.bam$", full.names = T)
	files_bed = list.files(path = Bedpath, pattern="\\.bed$", full.names = T)
	
	fbam = unlist(lapply(1:length(files_bam),function(i) sub('\\.bam$', '', basename(files_bam[i]))))
	fbed = unlist(lapply(1:length(files_bed),function(i) sub('\\.bed$', '', basename(files_bed[i]))))
	bambed = intersect(fbam,fbed)
	
    for (i in 1:length(bambed)){
        bamFile <- BamFile(paste0(Bampath,bambed[i],".bam"))
        bedfiles <- read.csv(paste0(Bedpath,bambed[i],".bed"), sep="\t", header=F)
        seqnames <- as.character(bedfiles[,1])
        starts <- ends <- bedfiles[,2]
        ranges=IRanges(starts,ends)
        w <- GRanges(seqnames,ranges)
        param <- ScanBamParam(which=w)
        res1 <- pileup(bamFile,scanBamParam=param)
        res = merge(res1,bedfiles, by=c(1,2), all.x=T)
        filename <- paste0(Outpath, bambed[i],".pileup")
		write.table(res, filename, row.names=F,sep='\t',quote =F)
    }
    
} else if (step==7){
	Outpath= args[2]
	library(rmarkdown)
	rmarkdown::render('TBias.Rmd')
 }
