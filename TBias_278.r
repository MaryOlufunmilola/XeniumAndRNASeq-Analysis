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
    
} else if (step == 7) {

    IOpath= args[2]
    
    library(reshape2)
    library(ggplot2)
	library(jsonlite)
    
	dict5 <- fromJSON("dict5.json", flatten=TRUE)
	dict6 <- fromJSON("dict6.json")	
	
	setwd(IOpath)
	for (i in 1:length(dict6)){
		a = data.frame(ref_alt = names(dict6[i][[1]]), tr = as.integer(sapply(dict6[i][[1]], '[[', 1)), ntr = as.integer(sapply(dict6[i][[1]], '[[', 2)))
		d6 = strsplit(dict5[i][[1]], "[\t]")
		xx = as.data.frame(do.call(rbind, d6))
		yy = xx[,3:5]
		colnames(yy) = c("ref_alt","tr","ntr")
		yy = yy[order(yy$ref_alt),]

		yy[,2] = as.numeric(levels(yy[,2]))[yy[,2]]
		yy[,3] = as.numeric(levels(yy[,3]))[yy[,3]]
		zz =melt(yy,id="ref_alt")
		x = split(zz, zz$ref_alt)
		
		if (any(sapply(x,nrow) < 4)){ 
			y = x[-which(sapply(x, nrow) < 4)]
		} else {
			y=x
		}

	    tests <- lapply(y, function(y) t.test(value ~ variable, y))
	    g = sapply(tests, "[[", "p.value")
	   
	    a <- a[which(a$ref_alt %in% names(g)),]
	    a$pval  <- lapply(g, round, 2)
	    a$Ref_alt = paste(a$ref_alt," ", "(p=",a$pval, ")")
	    
	    df <- melt(a[,-4])
	    x = sub('\\.pileup$', '', basename(names(dict5[i])))
	    png(file = paste0(x,'.png'))
	    print(ggplot(df, aes(factor(ref_alt), value, fill = variable)) + geom_bar(stat="identity", position = "dodge") +scale_fill_brewer(palette = "Set1") +
			   theme_bw() + labs(x="", y="", title=x) + theme(legend.position="top",axis.title.x=element_blank(), axis.text.x = element_text(angle=90)))
	    sink(paste0(x,'.csv'))
	    print(df)
	    sink()
	    dev.off()
	}
} else if (step==8){
	Outpath= args[2]
	library(rmarkdown)

	setwd(Outpath)
	rmarkdown::render('TBias.Rmd', output_dir =Outpath)
 }
