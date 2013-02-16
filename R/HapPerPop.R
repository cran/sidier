HapPerPop <-
function(readfile=T,sep=" ",header=F,inputFile=NA,input=NA,saveFile=T,Wname=NA,Iname=NA)
{
if (readfile==T)
{
input<-read.table(inputFile,sep=sep,header=header)
}
uhaplo<-sort(unique(input[,2]))
pops<-input[,1]
upops<-unique(pops)


matrizPresencia<-matrix(0,ncol=length(uhaplo),nrow=length(upops))
colnames(matrizPresencia)<-uhaplo
rownames(matrizPresencia)<-upops

for (i in 1:length(upops))
{
a<-input[which(upops[i]==input[,1]),2]
for(j in 1:length(a))
matrizPresencia[i,which(a[j]==colnames(matrizPresencia))]<-(1+matrizPresencia[i,which(a[j]==colnames(matrizPresencia))])
}

Pesos<-matrizPresencia
Interaccion<-matrizPresencia
Interaccion[which(matrizPresencia>0)]<-1

out<-list(c())
out[[1]]<-Pesos
out[[2]]<-Interaccion
names(out)<-c("Weighted","Interaction")

if(saveFile==T)
{
if(is.na(Wname))
Wname<-"Weighted.txt"
write.table(Pesos,file=Wname)

if(is.na(Iname))
Iname<-"Interaction.txt"
write.table(Interaccion,file=Iname)
}

out
}
