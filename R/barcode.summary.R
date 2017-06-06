barcode.summary<-function(dismat=NULL,save.distances=FALSE,folder.name="distance_matrices")
{

if(isSymmetric(dismat)==FALSE) stop("The input matrix must be symmetric")

mat<-dismat

sps<-row.names(mat)
kk<-sapply(sps,function(x){match(sps,x[1])})
SPS<-unique(sps)

Intra<-mat[which(kk==1 & lower.tri(mat))]
row.names(mat)[which(kk==1 & lower.tri(mat))]
colnames(mat)[which(kk==1 & lower.tri(mat))]
Inter<-mat[which(kk!=1 & lower.tri(mat))]
mat2<-mat
mat2[upper.tri(mat,diag=TRUE)]<-NA
for (i in 1:length(SPS))
  {
  SUB<-mat2[which(row.names(mat2)==SPS[1]),]
  if(i==1)
	{
	Intra<-list(sort(c(SUB[,which(colnames(SUB)==SPS[1])])))
	Inter<-list(sort(c(SUB[,which(colnames(SUB)!=SPS[1])])))
	}
  if(i>1)
	{
	Intra[[i]]<-sort(c(SUB[,which(colnames(SUB)==SPS[i])]))
	Inter[[i]]<-sort(c(SUB[,which(colnames(SUB)!=SPS[i])]))
	}
  }

outINTER<-t(sapply(Inter,summary))
row.names(outINTER)<-SPS
Ne<-sapply(Inter,length)

outINTRA<-t(sapply(Intra,summary))
row.names(outINTRA)<-SPS
Na<-sapply(Intra,length)

OUT<-list()
OUT[[1]]<-cbind(outINTRA,N=Na)
OUT[[2]]<-cbind(outINTER,N=Ne)

names(OUT)<-c("Intraspecific","Interspecific")
OUT
}
