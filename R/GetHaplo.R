GetHaplo <-
function(readfile=T,input=NA,align=NA,saveFile=T,outname="Haplotypes.txt",format="fasta",seqsNames=NA)
{
if(readfile==T)
align<-read.dna(file=input,format="fasta")
mat_alin<-as.matrix(as.character(align))
FH<-FindHaplo(readfile=readfile,input=input,align=align,saveFile=saveFile,outname=outname)
Huniques<-c()
U<-unique(FH[,2])
for(i in 1:length(U))
Huniques<-c(Huniques,which(FH[,2]==U[i])[[1]])
out<-align[Huniques,]

if(is.na(seqsNames[1])==F)
{
if(seqsNames[1]=="Inf.Hap")
dimnames(out)[[1]]<-U
if(seqsNames[1]!="Inf.Hap")
dimnames(out)[[1]]<-seqsNames
}

if(saveFile==T)
write.dna(out,file=outname,format=format)

print(paste(length(U)," different haplotypes found",ifelse(saveFile==T, {paste(", and saved in the file: \"",outname,"\"",sep="")},{paste(", but not saved in any file",sep="")}),sep=""),quote=F)


out

}
