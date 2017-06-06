inter.intra.plot<-function(dismat=NA, xlim=NULL,ylim=NULL, intra.col="gray",intra.density=0,intra.n=30,plot="N", inter.col="black",inter.density=0,inter.n=30,legend=TRUE,main="",xlab="Genetic distances",ylab=NULL)
{

if(isSymmetric(dismat)==FALSE) stop("The input matrix must be symmetric")

mat<-dismat
mat[upper.tri(mat,diag=T)]<-NA

intra<-list()
inter<-list()

RN<-matrix(rep(row.names(mat),ncol(mat)),byrow=F,ncol=ncol(mat))

CN<-matrix(rep(colnames(mat),nrow(mat)),byrow=T,ncol=ncol(mat))

INTRA<-sort(mat[which(RN==CN)])
INTER<-sort(mat[which(RN!=CN)])

aux1<-hist(INTER,plot=FALSE,n=inter.n)
aux2<-hist(INTRA,plot=FALSE,n=intra.n)

if(plot=="freq")
{
aux1$counts<-aux1$counts/sum(aux1$counts)
aux2$counts<-aux2$counts/sum(aux2$counts)
}

if(is.null(xlim)) xlim=c(min(range(aux1$breaks),range(aux2$breaks)),max(range(aux1$breaks),range(aux2$breaks)))

if(is.null(ylim))
ylim<-c(min(aux1$counts,aux2$counts),max(aux1$counts,aux2$counts))

if(is.null(ylab))
ifelse(plot=="freq", ylab<-"Frequency of pairwise comparisons",ylab<-"Number of pairwise comparisons")


plot(aux1,density=inter.density,xlim=xlim, ylim=ylim,main=main,xlab=xlab,ylab=ylab,border=inter.col)
plot(aux2,ylim=ylim,add=TRUE,xlim=xlim,density=intra.density, border=intra.col)
if(legend==TRUE)
legend("topright", legend=c("intra","inter"), pch=0,col=c(intra.col,inter.col))

out<-list()
out[[1]]<-INTRA
out[[2]]<-INTER
names(out)<-c("Intraspecific","Interspecific")

kk<-out
}

