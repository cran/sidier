single.network.module <-
function(dis,threshold=NA,refer2max=TRUE,out="module",save.file=FALSE,modFileName="Modules_summary.txt") # Meto cex.vertex y plot y get.coord
{
# require (igraph)
# require (network)

match.arg(out,c("module","network"))

if(is.na(threshold)==TRUE) print("ERROR: No threshold value defined")

dis<-as.matrix(dis)

## END na.rm

j<-threshold
dis2<-matrix(1,nrow=nrow(dis),ncol=ncol(dis))
row.names(dis2)<-row.names(dis)
ifelse(refer2max==TRUE,
lim<-max(dis)*j,
lim<-j)

fuera<-which(dis>lim)
dis2[fuera]<-0

G<-graph.adjacency(dis2)
A<-as.network.matrix(dis2)

comuni<-walktrap.community(G)
tab1<-matrix(nrow=nrow(dis2),ncol=2)
tab1<-as.data.frame(tab1)
tab1[,1]<-row.names(dis)
tab1[,2]<-comuni$membership
colores<-tab1[,2]
colnames(tab1)<-c("Node_label","Module")
if(save.file==TRUE)
write.table(file=modFileName,tab1,quote=FALSE,row.names=FALSE)

if(out=="module")
OUT<-tab1

if(out=="network")
OUT<-dis2
OUT
}

