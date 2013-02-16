FindHaplo <-
function(readfile=T,input=NA,align=NA,saveFile=T,outname="FindHaplo.txt")
{
library(ape)
if(readfile==T)
align<-read.dna(file=input,format="fasta")
mat_alin<-as.matrix(as.character(align))

seqs<-matrix(nrow=nrow(mat_alin),ncol=1)
for(i in 1:nrow(mat_alin))
seqs[i,]<-paste(mat_alin[i,],collapse="")

UH<-unique(seqs)
NH<-length(UH)
namesH<-matrix(nrow=NH,ncol=1)
for(i in 1:NH)
namesH[i,]<-(paste("H",i,sep=""))

out<-matrix(nrow=nrow(seqs),ncol=2)
colnames(out)<-c("Sequence.Name","Haplotype.Name")
out[,1]<-labels(align)

for(i in 1:NH)
out[which(seqs==UH[i,]),2]<-namesH[i,]

if(saveFile==T)
{
write.table(out,row.names=F,quote=F,file=outname)
}
out
}
