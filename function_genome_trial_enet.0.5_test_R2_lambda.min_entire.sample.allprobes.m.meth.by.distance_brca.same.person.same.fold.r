# 08/09/17
# cv R2 using all probes within 2Mbp from promoter region
# # genes 13646 
# elastic net alpha 0.5
library(GenomicRanges)
library(glmnet)
#source("divide.data.r")
#load("brca.data2.Rdata")
load("brca.data.non.missing.Rdata")
sample.id<-colnames(brca.data$meth)
patient.id <- substr(sample.id, 1, 12)

dup.tf <- duplicated(patient.id)|duplicated(patient.id, fromLast=TRUE)
dup.id <- patient.id[dup.tf]
dup.unique.id <- unique(dup.id)



load("probe.per.promoter_brca.Rdata")


promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

gene.lists <-promoter.range$gene.name[promoter.range$gene.name %in% rownames(brca.data$gene.exp)]

gene.list2 <- gene.lists[gene.lists %in% names(probe.per.promoter)]

n <- sum(ncol(brca.data$gene.exp))
curi<-1
nfolds<-10
lambda.rule <- "lambda.min"
#save(gene.list2, file="genome.list.Rdata")
#save(gene.list2, file="gene.list2.Rdata")
#gene.list2 <-gene.list2[1:20]
enhancer.range=2000000
library(parallel)
fn.R2.all.probe <- function(enhancer.range=2000000) {

  R2.all.probe.result <- mclapply(1:length(gene.list2),function(curi) {
   
    print(paste(curi,"th running",sep=""))
    promoter.range1 <- promoter.range[promoter.range$gene.name == gene.list2[curi]]
    
    candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),gene.name=promoter.range1$gene.name)
    probe.id <- unique(queryHits(findOverlaps(brca.data$meth.range, candi.enhancer.range1)))
    #probe.name.unique <- unique(c(promoter.probe.name, gene.body.probe.name))
    y=t(brca.data$gene.exp[rownames(brca.data$gene.exp) ==  gene.list2[curi],])
    x=t(brca.data$meth[ probe.id,])
   
   
   
    p<-ncol(data)-1
    # x1 = data[,-ncol(data)]
    set.seed(1234)
   
     
    ##############################################
    # enet.0.5
    cv.R2 <- rep(NA, 10)
    test.R2 <- rep(NA, 10)
    set.seed(1234)
    x <- apply(x, 2, as.double)
    for(i in 1:10) {
      #print(paste(i,"th running",sep=""))
       
      sam1<-  sample(c(0,1), size = length(dup.unique.id), replace=TRUE,prob=c(3/4,1/4))
      
      
      dup.unique.foldid<- data.frame(patient.id=dup.unique.id, foldid=sam1)
      sam2<- sample(c(0,1), size = sum(!dup.tf), replace=TRUE,prob=c(3/4,1/4))
      
      rest.foldid<- data.frame(patient.id=patient.id[!dup.tf], foldid=sam2)
      
      mergy.foldid <- rbind(dup.unique.foldid, rest.foldid)
      
      #foldids <- data.frame(sample.id=sample.id, mergy.foldid[match(patient.id, mergy.foldid$patient.id),])
      test.id <- mergy.foldid[match(patient.id, mergy.foldid$patient.id),"foldid"]
      
      
      x.train <- x[test.id==0,];x.test<-x[test.id==1,];y.train<-y[test.id==0];y.test<-y[test.id==1]
      cv.fit = cv.glmnet(x.train, y.train, alpha=0.5, keep=TRUE)
      id<-which(cv.fit$lambda == cv.fit$lambda.min)
      if(var(cv.fit$fit.preval[,id])!=0) {
        cv.R2[i] <- cor(cv.fit$fit.preval[,id],y.train)^2
      } else {
        cv.R2[i] <- 0
      }
      
      pred <- predict(cv.fit, newx=x.test, s="lambda.min", type = "response")
      if(var(pred)==0) {
        test.R2[i] <- 0
      } else {
        test.R2[i] <- cor(pred, y.test)^2
      }
    } 
    list(cv.R2=cv.R2, test.R2=test.R2)
  }, mc.cores=24)
  save(R2.all.probe.result, file=paste("all.probes.",c(enhancer.range/1000000),"M.expression.enet.0.5.cv.R2.m.meth.ngene=",length(gene.list2),".",lambda.rule,".entire.testR2.lung.Rdata",sep=""))
  
}




  
  

