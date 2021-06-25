get.majority.taxo<-function(filtered.taxo,verbose=TRUE){
x <- datos <- filtered.taxo
getProbTaxo.aux<-function(x,verbose=TRUE){
	X<-x[,(ncol(x)-6):ncol(x)]

	G1<-apply(X,2,function(x1){
		Tab1<-table(x1)
		noInfo<-which(names(Tab1)=="unidentified" | names(Tab1)=="Unidentified"  | names(Tab1)=="sp"  | names(Tab1)=="low_pident")

		ifelse(length(noInfo)!=0,{
		Tab2<-Tab1[-noInfo]
		Tab3<-sum (Tab1[noInfo])
		names(Tab3)<-"Uninformative"
		out<-c(sort(Tab2,decreasing=T),Tab3)},{out<-sort(Tab1,decreasing=T)})
		out})


		if(!is.na(match("Uninformative",lapply(G1,function(x) names(x[1])))))
		{
		keep<-(match("Uninformative",lapply(G1,function(x) names(x[1])))-1) 
		}else{keep<-7}

	# Si solo hay una taxo unica, la ponemos:
	#if(nrow(unique(X))==1)
	if(is.list(G1)==FALSE)
		{
		FINALtaxo<-X[1,]
		FINALtaxo[which(X[1,]=="unidentified" | X[1,]=="Unidentified"  | X[1,]=="sp"  | X[1,]=="low_pident")]<-"Uninformative"
	
		} else {

		while(!any(paste(lapply(G1,function(x) names(x[1]))[1:keep],collapse=";")==unique(apply(X[1:keep],1,paste,collapse=";"))))
		keep<-keep-1

		FINALtaxo<-unlist(lapply(G1,function(x) names(x[1])))[1:keep]
		if(length(FINALtaxo)<7)
			FINALtaxo<-c(FINALtaxo,rep("Uninformative",7-keep))
		}

	# Meter los numeros en otra columna:
	FINALtaxo<-c(FINALtaxo,matrix(paste(unlist(lapply(G1,paste,collapse="+")),collapse="|")))

	names(FINALtaxo)<-c("kingdom.final","phylum.final","class.final","order.final","family.final","genus.final","species.final","values")


	if(verbose==FALSE)
		{return(FINALtaxo)} else {
		kk<-list()
		kk[[1]]<-G1
		kk[[2]]<-FINALtaxo
		return(kk)
		}
	} # fin .aux	

## Aplicar .aux

seqs<-unique(datos$qseqid )
for (i in 1:length(seqs))
	{
	cat(paste(round(100*i/length(seqs),2),"%"),"\r")
	x<-datos[which(datos$qseqid==seqs[i]),]
	ifelse(i==1,
	OUT3<-t(matrix(c(unique(seqs[i]),getProbTaxo.aux(x,verbose=F)))),
	OUT3<-rbind(OUT3,c(unique(seqs[i]),getProbTaxo.aux(x,verbose=F))))
	#table(x[,12])
	}
colnames(OUT3)[1]<-"qseqid"
OUT3
}


