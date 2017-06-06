
rule<-function(summary=NULL,rule=NULL,stat.intra="max",stat.inter="min",pch.intra=16, pch.inter=16,pch.out=21,col.intra="gray",col.inter="black",col.out="black",label=F)
{

INTRA<-as.data.frame(summary$Intraspecific)
INTER<-as.data.frame(summary$Interspecific)

if(stat.intra=="max")
intraSTAT<-INTRA$Max.
if(stat.intra=="min")
intraSTAT<-INTRA$Min.
if(stat.intra=="median")
intraSTAT<-INTRA$Median.
if(stat.intra=="mean")
intraSTAT<-INTRA$Mean

if(stat.inter=="max")
interSTAT<-INTER$Max.
if(stat.inter=="min")
interSTAT<-INTER$Min.
if(stat.inter=="median")
interSTAT<-INTER$Median.
if(stat.inter=="mean")
interSTAT<-INTER$Mean

tabla<-data.frame(interSTAT,intraSTAT,Ratio=interSTAT/intraSTAT,(interSTAT/intraSTAT)>rule)

colnames(tabla)[c(1,2,4)]<-c(paste("Inter.",stat.intra,sep=""),paste("Intra.",stat.intra,sep=""),paste("Ratio>",rule,sep=""))
#	colnames(tabla)[4]<-paste(rule,"x",sep="")
row.names(tabla)<-row.names(INTER)


	
MediasX<-colMeans(tabla[which(tabla[,4]==TRUE),])[1:2]
PuntoMedioX<-(MediasX[2]-MediasX[1])/2+MediasX[1]

######
# Threshold to plot networks:
valoresX<-tabla[which(tabla[,4]==TRUE),][1:2]
MediasX<-colMeans(tabla[which(tabla[,4]==TRUE),])[1:2]
PuntoMedioX<-(MediasX[2]-MediasX[1])/2+MediasX[1]

# Graphic view of the gap:

plot(rep(1,2*nrow(tabla)),c(tabla[,1],tabla[,2]),pch=pch.out,col=col.out,xlab="",ylab="Genetic distance",xaxt="n")
 abline(h=PuntoMedioX)
 abline(h=MediasX[2],lty=2)
 abline(h=MediasX[1],lty=2)
points(rep(1,nrow(valoresX)),valoresX[,2],pch=pch.intra,col=col.intra)
points(rep(1,nrow(valoresX)),valoresX[,1],pch=pch.inter,col=col.inter)
legend("topright",pch=c(pch.inter,pch.intra,pch.out),col=c(col.inter,col.intra,col.out),legend=c(paste(rule,"x Inter. distances",sep=""),paste(rule,"x Intra. distances",sep=""),paste("Non-",rule,"x",sep="")),bg="white",box.col="white")
 box()

if(label==T)
if(any(tabla[,4]==TRUE))
{
text(y=tabla[which(tabla[,4]==TRUE),1],x=1,pos=4,label=row.names(tabla[which(tabla[,4]==TRUE),]),col=col.inter)
text(y=tabla[which(tabla[,4]==TRUE),2],x=1,pos=4,label=row.names(tabla[which(tabla[,4]==TRUE),]),col=col.intra)
}

text(x=0.7,y=MediasX[1]+0.01,labels=paste("Average",stat.inter,"inter"))
text(x=0.7,y=MediasX[2]+0.01,labels=paste("Average",stat.intra,"intra"))
text(x=0.7,y=PuntoMedioX+0.01,labels=paste(rule,"x Rule threshold",sep=""))

out<-matrix(nrow=3,ncol=1)
out[1]<-MediasX[2]
out[2]<-MediasX[1]
out[3]<-PuntoMedioX
out<-matrix(as.numeric(out))
colnames(out)<-"Genetic distance"
row.names(out)<-c(paste("Average",stat.intra,"intra"),paste("Average",stat.inter,"inter"),paste(rule,"x Rule threshold",sep=""))

percent<-round(length(which(tabla[,4]==TRUE))/nrow(tabla)*100,2)
cat(paste("\nThreshold estimated using a subset of ",length(which(tabla[,4]==TRUE))," species, representing ",percent,"% of the studied species\n",sep=""))

out<-list(out,tabla)
names(out)<-c("Rule","Summary")
out
}



