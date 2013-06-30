perc.thr <-
function(dis,range=seq(0,1,0.01),ptPDF=TRUE,ptPDFname="PercolatedNetwork.pdf",estimPDF=TRUE,estimPDFname="PercThr Estimation.pdf",estimOutfile=TRUE, estimOutName="PercThresholdEstimation.txt",appendOutfile=TRUE,plotALL=FALSE,bgcol="white",label.col="black",label=colnames(dis),modules=FALSE,moduleCol=NA,modFileName="Modules_summary.txt")
{
require(igraph)
require(network)

salida<-matrix(nrow=length(range),ncol=3)
colnames(salida)<-c("Threshold","<s>","Clusters")

for (j in range)
{
print(paste("Threshold value:",j,"  Range to test: from ",min(range)," to ",max(range),sep=""))
dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
lim<-max(dis)*j
fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

Res<-clusters(G)
noGrande<-sort(Res$csize)[-length(sort(Res$csize))]
N<-sum(noGrande)
repes<-unique(noGrande[which(duplicated(noGrande))])
if (Res$no>1)
{
if (length(repes)>0)
{
n<-c()
for (i in 1:length(repes))
n<-c(n,length(which(noGrande==repes[i])))
sum1<-repes^2*n

noUnic<-c()
for (i in 1:length(repes))
noUnic<-c(noUnic,which(noGrande==repes[i]))
unicos<-noGrande[-noUnic]
sum2<-unicos^2

SUM<-sum(sum1)+sum(sum2)
S<-SUM/N
}

if (length(repes)==0)
S<-sum(noGrande^2)/N
}
if (Res$no==1)
S<-1
#print(c(j,S,Res$no))

salida[which(range==j),1]<-j
salida[which(range==j),2]<-S
salida[which(range==j),3]<-Res$no

if(is.null(colnames(dis)))
label<-c(1:ncol(dis))

if(plotALL==TRUE)
	{
pdf(file=paste("Threshold=",j,".pdf",sep=""))
		if(modules==T)
		{
		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
		colores<-tab1[,2]
		bgcol<-colores
		}

	plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5,interactive=F,label.pos=5,label.col=label.col,label.cex=0.8,main=paste("Threshold=",j,sep=" "))
dev.off()
#	dev.copy2pdf(file=paste("Threshold=",j,".pdf",sep=""))
	}
}

print("Preparing outfiles")

sal<-salida[,-3]
if(estimPDF==TRUE)
{
pdf(file=paste("Threshold=",j,".pdf",sep=""))
plot(sal,type="l")
points(sal)
#dev.copy2pdf(file=estimPDFname)
dev.off()
}

if(estimPDF==FALSE)
{
plot(sal,type="l")
points(sal)
#dev.copy2pdf(file=estimPDFname)
}

if(estimOutfile==TRUE)
write.table(salida,file=estimOutName,append=appendOutfile,row.names=FALSE)

out<-list(c())
out[[1]]<-salida
out[[2]]<-salida[(max(which(salida[,2]>1))+1),1]
names(out)<-c("Summary","Estimated Percolation Threshold")

#
j<-salida[(max(which(salida[,2]>1))+1),1]
dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
row.names(dis2)<-row.names(dis)
lim<-max(dis)*j
fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

		if(modules==T)
		{
		comuni<-walktrap.community(G)
		tab1<-matrix(nrow=nrow(dis2),ncol=2)
		tab1<-as.data.frame(tab1)
		tab1[,1]<-label
		tab1[,2]<-comuni$membership
		colo<-colors()[sample(c(1,23,25:152,203:259,361:657),length(unique(tab1[,2])))]
		if(is.character(moduleCol[1])==T)
		colo<-moduleCol

		tab1[which(tab1[,2]==1),3]<-colo[1]
		if(length(unique(tab1[,2]))>1)
		for(i in 2:length(unique(tab1[,2])))
		tab1[which(tab1[,2]==i),3]<-colo[i]
		colnames(tab1)<-c("Node_label","Module","Node_colour")

		bgcol<-tab1[,3]
		out[[3]]<-tab1
		names(out)<-c("Summary","Estimated Percolation Threshold","Module")
		
		write.table(file=modFileName,tab1,quote=F,row.names=FALSE)
		}

plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5,interactive=F,label.pos=5,label.col=label.col,label.cex=0.8,main=paste("Threshold=",j,sep=" "))

if(ptPDF==TRUE)
{
pdf(file=ptPDFname)
plot.network(A,vertex.col=as.matrix(bgcol),label=label,usearrows=0,vertex.cex=2.5,interactive=F,label.pos=5,label.col=label.col,label.cex=0.8,main=paste("Threshold=",j,sep=" "))
dev.off()
#dev.copy2pdf(file=ptPDFname)
}

out
}
