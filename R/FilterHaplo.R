FilterHaplo<-function(inputFile=NA,align=NA,Nmin=0,Nmax=NULL,saveFile=FALSE,outname="FilterHaplo.txt"){

if(is.na(inputFile)==TRUE&is.na(align[1])==TRUE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==FALSE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==TRUE)
align<-read.dna(file=inputFile,format="fasta")

FH<-FindHaplo(align=align,saveFile=FALSE)
HaploFreq<-table(FH[,2])

if(is.null(Nmax)) Nmax=max(HaploFreq)+100
new.align<- align
new.FH<-FindHaplo(align=new.align,saveFile=FALSE)
new.HaploFreq<-table(new.FH[,2])
old.names<-names(HaploFreq)
names(new.HaploFreq)<-old.names

if(length(which(HaploFreq<Nmin))!=0 & length(which(HaploFreq>Nmax))==0)
 {
 haplos.lower<-names(which(HaploFreq<Nmin))
remove.hap<-c()
 for (i in haplos.lower)
	remove.hap<-c(remove.hap,which(FH[,2]==i))

new.align<- align[-c(remove.hap),]
new.FH<-FindHaplo(align=new.align,saveFile=FALSE)
new.HaploFreq<-table(new.FH[,2])
old.names<-names(HaploFreq)[-c(match(c(haplos.lower), names(HaploFreq)))]
names(new.HaploFreq)<-old.names
 }

if(length(which(HaploFreq<Nmin))==0 & length(which(HaploFreq>Nmax))!=0)
 {
 haplos.higher<-names(which(HaploFreq>Nmax))
 if(exists("remove.hap")==FALSE)
	 remove.hap<-c()
 for (i in haplos.higher)
	remove.hap<-c(remove.hap,which(FH[,2]==i))

new.align<- align[-c(remove.hap),]
new.FH<-FindHaplo(align=new.align,saveFile=FALSE)
new.HaploFreq<-table(new.FH[,2])
old.names<-names(HaploFreq)[-c(match(c(haplos.higher), names(HaploFreq)))]
names(new.HaploFreq)<-old.names
 }

if(length(which(HaploFreq<Nmin))!=0 & length(which(HaploFreq>Nmax))!=0)
 {
 haplos.lower<-names(which(HaploFreq<Nmin))
remove.hap<-c()
 for (i in haplos.lower)
	remove.hap<-c(remove.hap,which(FH[,2]==i))

 haplos.higher<-names(which(HaploFreq>Nmax))
 if(exists("remove.hap")==FALSE)
	 remove.hap<-c()
 for (i in haplos.higher)
	remove.hap<-c(remove.hap,which(FH[,2]==i))

new.align<- align[-c(remove.hap),]
new.FH<-FindHaplo(align=new.align,saveFile=FALSE)
new.HaploFreq<-table(new.FH[,2])
old.names<-names(HaploFreq)[-c(match(c(haplos.lower,haplos.higher), names(HaploFreq)))]
names(new.HaploFreq)<-old.names
 }



cat("\n Haplotype frequency in the original alignment\n")
print(HaploFreq)
cat("\n Haplotype frequency in the new alignment\n\n")
print(new.HaploFreq)
new.align
}

