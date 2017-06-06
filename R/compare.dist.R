compare.dist<-function(distr1=NULL,distr2=NULL,N=50,normalize=TRUE,main=NA,col1="gray",col2="black", col.border1="gray",col.border2="black",col.line1="gray",col.line2="black",Ylab=c("Abundance","Abundance","Success"), Xlab=c("data1","data2","Threshold"))
{
data1<- na.omit(distr1)
data2<- na.omit(distr2)

	cut.lines<-function(A,B,C,D,plot=FALSE) # A and B are points for line 1 and C-D are for line 2
	{
	Ax<-A[1]
	Bx<-B[1]
	Cx<-C[1]
	Dx<-D[1]
	Ay<-A[2]
	By<-B[2]
	Cy<-C[2]
	Dy<-D[2]
	mAB=(By-Ay)/(Bx-Ax)
	nAB=Ay-mAB*Ax

	mCD=(Dy-Cy)/(Dx-Cx)
	nCD=Cy-mCD*Cx

	X=(nCD-nAB)/(mAB-mCD)
	Y= mAB*X+nAB

	out<-matrix(c(X,Y),ncol=2,dimnames=list(c(),c("X","Y")))


	if (plot==TRUE)
		{
		plot(c(Ax,Bx,Cx,Dx),c(Ay,By,Cy,Dy),xlim=c(min(Ax,Bx)-5,max(Ax,Bx)+5))
		abline(a=nAB,b=mAB)
		abline(a=nCD,b=mCD)
		points(X,Y,pch=16,col="red")
		}
	out
	}

#ifelse(min(distr1)<min(distr2),{data1<-distr1;data2<-distr2},ifelse(min(distr1)==min(distr2),{data1<-distr1;data2<-distr2;warning("try changing distr1 and distr2")},{data1<-distr1;data2<-distr2}))


MIN<-min(data1)
MAX<-max(data2)
STEP<-(MAX-MIN)/(N-3)
secu<-seq(MIN-STEP,MAX+STEP,STEP)

abo<-c()
for (i in 1:length(secu))
abo<-c(abo,length(which(data1<secu[i])))

bel<-c()
for (i in 1:length(secu))
bel<-c(bel,length(which(data2>secu[i])))

if(normalize==TRUE)
	{
	abo<-100*abo/max(abo)
	bel<-100*bel/max(bel)
#	Ylab[3]<-paste(Ylab[3],"(%)")
	}

layout(matrix(c(1:3),nrow=3))
if(is.na(main[1]))
{
hist(data1,n=N,xlim=c(MIN,MAX),col=col1,border=col.border1,xlab=Xlab[1],ylab=Ylab[1])
hist(data2,n=N,xlim=c(MIN,MAX),col=col2,border=col.border2,xlab=Xlab[2],ylab=Ylab[2])
}
if(is.na(main[2])==FALSE)
{
hist(data1,n=N,xlim=c(MIN,MAX),col=col1,border=col.border1,main=main[1],xlab=Xlab[1],ylab=Ylab[1])
hist(data2,n=N,xlim=c(MIN,MAX),col=col2,border=col.border2,main=main[2],xlab=Xlab[2],ylab=Ylab[2])
}

plot(x=secu,y=abo,type="l",xlim=c(MIN,MAX),col=col.line1,ylab=Ylab[3],xlab=Xlab[3])
points(x=secu,y=bel,type="l",xlim=c(MIN,MAX),col=col.line2)

both<-cbind(abo,bel)

iguales<-apply(both,1,function(x){x[1]==x[2]})
if(sum(iguales)==1) {out1<-secu[which(iguales==TRUE)]; out2<-both[which(iguales==TRUE),1]; out<-matrix(c(out1,out2),ncol=2,dimnames=list(c(),c("X","Y")))}

# Si hay muchos iguales, cogemos el punto medio de ellos
#if(sum(iguales)>1) {out1<-(secu[max(which(iguales==TRUE))]-secu[min(which(iguales==TRUE))])/2;out2<-both[max(which(iguales==TRUE))]-secu[min(which(iguales==TRUE)),1]; out<-matrix(c(out1,out2),ncol=2,dimnames=list(c(),c("X","Y")))}

if(sum(iguales)>1) {out1<-mean(c(min(secu[which(iguales==TRUE)]),max(secu[which(iguales==TRUE)]))); out2<-both[which(secu==out1),1]; out<-matrix(c(out1,out2),ncol=2,dimnames=list(c(),c("X","Y")))}

if(sum(iguales)==0)
	{
	aux1<-max(which(apply(both,1,function(x){x[1]<x[2]})==TRUE))
	aux2<-min(which(apply(both,1,function(x){x[1]>x[2]})==TRUE))

	A=c(secu[aux1],abo[aux1])
	B=c(secu[aux2],abo[aux2])
	C=c(secu[aux1],bel[aux1])
	D=c(secu[aux2],bel[aux2])
	out<-cut.lines(A,B,C,D,plot=FALSE)
	}

points(c(MIN,out[1]),c(out[2],out[2]),type="l",lty=2)
points(c(out[1],out[1]),c(out[2],0),type="l",lty=2)

text(x=MIN,y=out[2],label=round(out[2],2),pos=1)
text(x=out[1],y=0,label=round(out[1],2),pos=2)

colnames(out)<-c("Threshold","Success(%)")
out

}

