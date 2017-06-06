barcode.quality<-function(dismat=NA,threshold=NA,refer2max=FALSE,save.file=FALSE,modFileName="Modules_summary.txt",verbose=FALSE,output="list")
{

T0<-single.network.module(dismat,threshold,refer2max,out="module",save.file=save.file,modFileName=modFileName)

success<-function(x)
{
MEDIAS<-c()
out<-c()

	for (i in 1:nrow(x))
	{
	this.group<-x[which(x[,2]==x[i,2]),]
	OK<-length(which(this.group[,1]==x[i,1]))
	NOTlink<-length(which(this.group[,1]!=x[i,1]))
	OKall<-nrow(x[which(x[,1]==x[i,1]),])

	out<-c(out,OK/(OKall+NOTlink))
	}
	MEDIAS<-c(MEDIAS,mean(out))

matrix(MEDIAS,dimnames=list("Qvalue",""))
}

acc.prec<-function(T0)
	{
	x<-unique(T0[,2])
	accu<-c()
	prec<-c()

	TRUE_ok<-0
	TRUE_no<-0
	FALSE_ok<-0
	FALSE_no<-0
	
		for (i in x)
		{
		this.group<-T0[which(T0[,2]==i),]
		this.sp<- names(sort(table(this.group[,1]),decreasing=TRUE)[1])
		TRUE_ok<-TRUE_ok+sort(table(this.group[,1]),decreasing=TRUE)[1]
		TRUE_no<-c(TRUE_no+length(which(T0[,1]!=this.sp & T0[,2]!=i)))
		FALSE_no<-c(FALSE_no+length(which(T0[,1]==this.sp & T0[,2]!=i)))
		FALSE_ok<-c(FALSE_ok+sum(sort(table(this.group[,1]),decreasing=TRUE)[-1]))
		}
	
	accu<-(TRUE_ok+TRUE_no)/(TRUE_ok+TRUE_no+FALSE_no+FALSE_ok)
	prec<-TRUE_ok/(TRUE_ok+FALSE_ok)
	Fscore<-TRUE_ok/(TRUE_ok+FALSE_no+FALSE_ok)

matrix(c(accu,prec,Fscore),nrow=3,dimnames=list(c("Accuracy","Precision","Fscore"),""))
}

OUT<-rbind(success(T0),acc.prec(T0))

if(verbose==TRUE) 
	{
	N.modules<-max(T0[,2])
	N.sp.mod<-c()
	Modus<-unique(T0[,2])
		for (i in Modus)
		N.sp.mod<-c(N.sp.mod,length(unique(T0[which(T0[,2]==i),1])))
	N.sp.mod.MAX<-max(N.sp.mod)
	N.sp.mod.MED<-mean(N.sp.mod)
	N.sp.mod.1<-length(which(N.sp.mod==1))
	
	N.species<-length(unique(T0[,1]))
	N.mod.sp<-c()
	Spes<-unique(T0[,1])
		for (i in Spes)
		N.mod.sp<-c(N.mod.sp,length(unique(T0[which(T0[,1]==i),2])))
	N.mod.sp.MAX<-max(N.mod.sp)
	N.mod.sp.MED<-mean(N.mod.sp)
	N.mod.sp.1<-length(which(N.mod.sp==1))

## metido nuevo
	kk4<-unique(T0[unique(order(T0[,2])),])
	kk5<-kk4[match( which(table(kk4[,2])==1), kk4[,2]),]

nomb<-unique(kk5[,1])
 unos<-c()
 for (i in nomb)
 unos<-c(unos,length(which(kk5[,1]==i)))
 N.mod.FIT.sp<-length(which(unos==1))
###




	out<-list()
	out[[1]]<-N.modules
	out[[2]]<-rbind(N.sp.mod.1,N.sp.mod.MAX,N.sp.mod.MED)
	out[[3]]<-N.species
	out[[4]]<-rbind(N.mod.sp.1,N.mod.sp.MAX,N.mod.sp.MED)
	out[[5]]<-N.mod.FIT.sp
	out[[6]]<-OUT

	names(out)<-c("Number.of.modules","Number.of.species.per.module","Number.of.species", "Number.of.modules.per.species","Number.of.modules.fitting.defined.species","Quality.estimates")
	OUT<-out

	if(output=="matrix")
		{
		OUT<-rbind(N.modules,N.sp.mod.1,N.sp.mod.MAX,N.sp.mod.MED,N.species,N.mod.sp.1,N.mod.sp.MAX,N.mod.sp.MED,N.mod.FIT.sp,success(T0),acc.prec(T0))
#		row.names(OUT)<-c("Number.of.modules", "Max.Number.of.species.per.module" ,"Mean.Number.of.species.per.module", "Number.of.modules.with.1.species","Number.of.species", "Max.Number.of.modules.per.species", "Mean.Number.of.modules.per.species", "Number.of.species.in.1.module","Q-value","Accuracy","Precision","Fscore")
		}
	}

OUT
}


