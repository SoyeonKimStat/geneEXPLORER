test<- read.table("HumanMethylation450",header=TRUE, row.names=1,check.names=FALSE)
m <- function(beta) {
  log2(beta/(1-beta))
}

meth.brca.m <- apply(test,2, m)

#save(meth.brca.m, file="meth.brca.m.Rdata")

#load("meth.brca.m.Rdata")
any(is.infinite(meth.brca.m))
#meth.lugn.m.org <- meth.brca.m
#meth.brca.m[is.infinite(meth.brca.m.org)] <- NA
#meth.brca.m["cg18106189",]
any(is.na(meth.brca.m))
any.na<- apply(meth.brca.m, 1, function(x) any(is.na(x)))

sum(any.na)/nrow(meth.brca.m) #25.1% missing or infinite value

n <- ncol(meth.brca.m)
percent.na<- apply(meth.brca.m, 1, function(x) sum(is.na(x))/n)
hist(percent.na[percent.na!=0])

with.na<- percent.na[percent.na!=0]
sum(with.na <=0.3)
sum(with.na <=0.2)
impute.ind <- percent.na!=0&percent.na<=0.2
which.na<-apply(meth.brca.m[impute.ind,], 1, function(x) which(is.na(x)))
which.na2<-is.na(meth.brca.m[impute.ind,])
sum(impute.ind)
library(REMP)
temp <- grooMethy(meth.brca.m[impute.ind,])
imputed<-getM(temp)


hist(imputed[which.na2])
hist(imputed[!which.na2])

any(is.na(imputed))

meth.brca.m.imputed <- meth.brca.m
meth.brca.m.imputed[impute.ind,] <- imputed

any.na<- apply(meth.brca.m.imputed, 1, function(x) any(is.na(x)))

sum(any.na)/nrow(meth.brca.m.imputed) #18.5% missing

#save(meth.brca.m.imputed, file="meth.brca.m.imputed.Rdata")

meth.brca.m.imputed.no.na <- meth.brca.m.imputed[!any.na,] # no missing value only
#save(meth.brca.m.imputed.no.na, file="meth.brca.m.imputed.no.na.Rdata")

##################################################
# gene expression
gene.exp0<- read.table("HiSeqV2",header=TRUE, row.names=1,check.names=FALSE)

gene.exp.mean <- rowMeans(gene.exp0)
sum(gene.exp.mean<=1)
high.gene.exp <- gene.exp0[gene.exp.mean>1,]
meth0<-meth.brca.m.imputed.no.na[,colnames(meth.brca.m.imputed.no.na) %in% colnames(high.gene.exp)]
gene.exp <- high.gene.exp[,colnames(meth0)]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, 
                       what = "Locations")
probe.range <- GRanges(seqnames = probe$chr, ranges = IRanges(probe$pos, 
                                                              width = 1, names = rownames(probe)), strand = probe$strand, 
                       name = rownames(probe))

id <- match(rownames(meth0), names(probe.range))
sort.non.na.id <- sort(id[!is.na(id)])
meth.range <- probe.range[sort.non.na.id]

meth<- meth0[names(meth.range),]

brca.data <- list(meth=meth, gene.exp=gene.exp, meth.range=meth.range)
#save(brca.data, file="brca.data2.Rdata" )
######################################
# promoter.per.probe
promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

#load("brca.data2.Rdata")

probe.range <- brca.data$meth.range
ids <- unique(subjectHits(findOverlaps(probe.range, promoter.range)))
sids <- sort(ids)
library(parallel)
probe.per.promoter<-mclapply(1:length(sids), function(curi){
  num<-sids[curi]
  probe.range[queryHits(findOverlaps(probe.range, 
                                     promoter.range[num]))]
},mc.cores=20)
names(probe.per.promoter) <- promoter$gene.name[sids]
save(probe.per.promoter, file="probe.per.promoter_brca.Rdata")
##############################################################
#synonym search
#load("brca.data2.Rdata")

promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)


source("synonym.search.r")
syn.name<- synonym.search(promoter.range$gene.name,rownames(brca.data$gene.exp))

org.name<-rownames(brca.data$gene.exp)


brca.data$gene.exp<- as.matrix(brca.data$gene.exp)
rownames(brca.data$gene.exp) <- syn.name
sum(!org.name %in% syn.name) #1492
save(brca.data,file="brca.data.non.missing.Rdata")