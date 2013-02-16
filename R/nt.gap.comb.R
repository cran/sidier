nt.gap.comb <-
function(DISTnuc=NA,DISTgap=NA,range=seq(0,1,0.1),method="Corrected",saveFile=TRUE)
{
library(ape)

CORnuc<-DISTnuc/(max(DISTnuc))
CORgap<-DISTgap/(max(DISTgap))

OUTuncor<-list(c())
OUTcor<-list(c())

if(method=="Uncorrected"|method=="Both")
{
for(ite1 in 1:length(range))
{
alfa<-range[ite1]
DIST<-(1-alfa)*DISTnuc+alfa*DISTgap
if(saveFile==T)
write.table(DIST,file=paste("DistanceMatrixUncorrectedAlfa",alfa,sep="_"))
OUTuncor[[ite1]]<-DIST
}
}

if(method=="Corrected"|method=="Both")
{
for(ite1 in 1:length(range))
{
alfa<-range[ite1]
COR<-(1-alfa)*CORnuc+alfa*CORgap
if(saveFile==T)
write.table(COR,file=paste("DistanceMatrixCorregidaAlfa",alfa,sep="_"))
OUTcor[[ite1]]<-COR
}
}

if(method=="Uncorrected"|method=="Both")
names(OUTuncor)<-paste("Alpha=",range)

if(method=="Corrected"|method=="Both")
names(OUTcor)<-paste("Alpha=",range)

OUT<-list(OUTuncor,OUTcor)
names(OUT)<-c("Uncorrected","Corrected")

OUT

}
