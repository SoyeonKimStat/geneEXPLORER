# 01/23/21
# install GenomicRanges, glmnet, data.table, R.utils, and REMP
# input study abbreviation
# default range is +-10Mb from promoter regions (enhancer.range)

library(GenomicRanges)
library(glmnet)
#source("divide.data.r")
study <- "BRCA"

file <- paste("https://tcga.xenahubs.net/download/TCGA.",study,".sampleMap/HumanMethylation450.gz",sep="")
download.file(file, paste("HumanMethylation450_",study,".gz", sep=""))
library(R.utils)
gunzip(paste("HumanMethylation450_",study,".gz", sep=""))
file2 <- paste("https://tcga.xenahubs.net/download/TCGA.",study,".sampleMap/HiSeqV2.gz",sep="")

download.file(file2, paste("HiSeqV2_",study,".gz", sep=""))
gunzip(paste("HiSeqV2_",study,".gz", sep=""))

source("read.data.preprocessing_function_infinite.r")
data.processing(study)

load(paste(study,".data.Rdata", sep=""))
load(paste("probe.per.promoter_",study,".Rdata", sep=""))
sample.id<-colnames(data$meth)
patient.id <- substr(sample.id, 1, 12)

dup.tf <- duplicated(patient.id)|duplicated(patient.id, fromLast=TRUE)
dup.id <- patient.id[dup.tf]
dup.unique.id <- unique(dup.id)
#test<-foldids[patient.id%in%dup.unique.id ,]
#test[order(test$patient.id),]


promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

gene.lists <-promoter.range$gene.name[promoter.range$gene.name %in% rownames(data$gene.exp)]

gene.list2 <- gene.lists[gene.lists %in% names(probe.per.promoter)]

n <- sum(ncol(data$gene.exp))
curi<-1
nfolds<-10
lambda.rule <- "lambda.min"
#save(gene.list2, file="genome.list.Rdata")
#save(gene.list2, file="gene.list2.Rdata")
#gene.list2 <-gene.list2[1:20]
enhancer.range=1e7
library(parallel)
fn.cv.R2 <- function(seq.num, k) {

  sub.gene.list <- gene.list2[seq.num]
  cv.R2.all.probe.result <- mclapply(1:length(seq.num),function(curi) {
    
     
    print(paste(curi,"th running",sep=""))
    promoter.range1 <- promoter.range[promoter.range$gene.name == sub.gene.list[curi]]
    
    candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),gene.name=promoter.range1$gene.name)
    probe.id <- unique(queryHits(findOverlaps(data$meth.range, candi.enhancer.range1)))
    #probe.name.unique <- unique(c(promoter.probe.name, gene.body.probe.name))
    y=t(data$gene.exp[rownames(data$gene.exp) ==  sub.gene.list[curi],,drop=FALSE])
    x=t(data$meth[ probe.id,])
    
    
    
    p<-ncol(data)-1
    # x1 = data[,-ncol(data)]
    set.seed(1234)
    sam1<- sample(rep(seq(nfolds), length = length(dup.unique.id)))
    
    dup.unique.foldid<- data.frame(patient.id=dup.unique.id, foldid=sam1)
    sam2<- sample(rep(nfolds:1, length = sum(!dup.tf)))
    
    rest.foldid<- data.frame(patient.id=patient.id[!dup.tf], foldid=sam2)
    
    mergy.foldid <- rbind(dup.unique.foldid, rest.foldid)
    
    #foldids <- data.frame(sample.id=sample.id, mergy.foldid[match(patient.id, mergy.foldid$patient.id),])
    foldid <- mergy.foldid[match(patient.id, mergy.foldid$patient.id),"foldid"]
    
    
    x=apply(x,2,as.double)
    ny <- ncol(y)
    if(ny ==1) {
      y <- as.numeric(y)
      cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5, foldid=foldid)
      #coef(cv.fit, s=lambda.rule)
      id<-which(cv.fit$lambda == cv.fit$lambda.min)
      if(var(cv.fit$fit.preval[,id])!=0) {
        R2 <- cor(cv.fit$fit.preval[,id],y)^2
      } else {
        R2 <- 0
      }
    } else {
      
      sR2 <- rep(0, ny)
      for(i in 1:ny) {
        y.m <- as.numeric(y[,i])
        cv.fit <- cv.glmnet(x, y.m, keep=TRUE, alpha=0.5, foldid=foldid)
        id<-which(cv.fit$lambda == cv.fit$lambda.min)
        
        if(var(cv.fit$fit.preval[,id])!=0) {
          sR2[i] <- cor(cv.fit$fit.preval[,id],y.m)^2
        } else {
          sR2[i] <- 0
        }
        
      } 
      R2 <- max(sR2)
    }
    names(R2) <- sub.gene.list[curi]
    list(cv.fit=cv.fit, R2=R2)
    }, mc.cores=50)
  #save(test.R2, file=paste("test.R2.GSE39004.m.meth.ngene=",length(test.R2),".testR2.Rdata",sep=""))
  save(cv.R2.all.probe.result, file=paste("cv.R2.",study,".",c(enhancer.range/1000000),"M.expression.enet.0.5.m.meth.",k,"th.running.Rdata",sep=""))
}

all.cv.R2 <- rep(NA,length(gene.list2))
names(all.cv.R2) <- gene.list2
train.models <- vector("list", length=length(gene.list2))
names(train.models) <- gene.list2
seq.num <- 1: length(gene.list2)

k<-1
while(length(seq.num) !=0) {
  print(paste(k,"th.running",sep=""))
  
  fn.cv.R2(seq.num, k)
  
  load(paste("cv.R2.",study,".",c(enhancer.range/1000000),"M.expression.enet.0.5.m.meth.",k,"th.running.Rdata",sep=""))
  R2 <- unlist(lapply(cv.R2.all.probe.result, function(x) x$R2))
  cv.fits <-lapply(cv.R2.all.probe.result, function(x) x[1])
   
  all.cv.R2[match(names(R2), names(all.cv.R2))] <- R2
  train.models[match(names(R2), names(all.cv.R2))] <- cv.fits
  k <- k+1
  seq.num <- which(is.na(all.cv.R2))
  
}

all.cv.R2 <- all.cv.R2[!is.na(all.cv.R2)]
#pred.gene.exp <- pred.gene.exp[!is.na(all.cv.R2), ]
n <- ncol(data$meth)
save(list=c('all.cv.R2', 'train.models'), file=paste("cv.R2.",study,".",c(enhancer.range/1000000),"M.expression.enet.0.5.m.meth.train.models.Rdata",sep="")) 

#fit <- train.models[[which(names(train.models) == "CHMP1B")]]

#coef(fit$cv.fit)