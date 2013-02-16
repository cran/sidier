pop.dist <-
function(DistFile=T,inputDist=NA,distances=NA,HaploFile=T,inputHaplo=NA,Haplos=NA,outType="O",logfile=TRUE,saveFile=TRUE,NameIni=NA,NameEnd=NA)
{
ifelse(HaploFile==T,
Haplos<-read.table(file=inputHaplo),
Haplos<-data.frame(Haplos))

pops<-unique(row.names(Haplos))

if(DistFile==T)
distances<-read.table(inputDist)
#row.names(distances)<-substr(row.names(distances),1,3)
#colnames(distances)<-substr(colnames(distances),1,3)


POPdist<-matrix("NA",nrow=length(pops),ncol=length(pops))
row.names(POPdist)<-pops
colnames(POPdist)<-pops

for (k in 1:nrow(POPdist))
POPdist[k,k]<-0

for(k in 1:(nrow(Haplos)-1))
for(l in (k+1):nrow(Haplos))
{
pop1<-Haplos[k,]
pop2<-Haplos[l,]

Hpop1<-rep(colnames(pop1)[which(pop1!=0)],pop1[which(pop1!=0)])
Hpop2<-rep(colnames(pop2)[which(pop2!=0)],pop2[which(pop2!=0)])

DIST1_2<-c()
for(i in 1:length(Hpop1))
for(j in 1:length(Hpop2))
DIST1_2<-c(DIST1_2,distances[which(row.names(distances)==substr(Hpop1[i],NameIni,NameEnd)),which(colnames(distances)==substr(Hpop2[i],NameIni,NameEnd))])
if(length(DIST1_2)==0) DIST1_2<-0
if(outType=="7"|outType=="O")
POPdist[k,l]<-mean(DIST1_2)
if(outType=="L"|outType=="O")
POPdist[l,k]<-mean(DIST1_2)

row.names(POPdist)<-pops
colnames(POPdist)<-pops
}
if(logfile==T)
write.table(c(paste("Among haplotypes distance matrix used:",inputDist),paste("Haplotypes per population matix used",inputHaplo)),file=paste(inputDist,"_PopopulationDistances.r.txt.log",sep=""))
if(saveFile==T)
write.table(POPdist,file=paste(inputDist,"_PopulationDistances.r.txt",sep=""),na="")
#matrix(as.numeric(POPdist),nrow=nrow(POPdist))
#POPdist
as.data.frame(POPdist,nrow=nrow(POPdist))
}
