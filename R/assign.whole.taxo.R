assign.whole.taxo<-function(BLAST){
BEST<-BLAST
get.sp<-function(x)
	{
	titulosUNITE2<-sapply(x,function(x){strsplit(x,"|",fixed=T)[[1]][length(strsplit(x,"|",fixed=T)[[1]])]})
	titulosUNITE4<-sapply(titulosUNITE2,function(x){strsplit(x,";")})
	GETsp<-	paste(unlist(strsplit(substr(titulosUNITE4[[1]][7],4,nchar(x)),"_")),collapse=" ")
	GETwhole<-t(as.matrix(c(substr(titulosUNITE4[[1]][1:6],4,1000),strsplit(substr(titulosUNITE4[[1]][7],4,nchar(x)),"_")[[1]][2])))
	
	as.matrix(cbind(GETsp,GETwhole))
	}

OUT<-cbind(BEST,"Organism"=NA,"kingdom.whole"=NA,"phylum.whole"=NA,"class.whole"=NA, "order.whole"=NA,"family.whole"=NA,"genus.whole"=NA,"species.whole"=NA)

for (i in 1:nrow(BEST))
	{
	cat(paste(round(100*i/nrow(BEST),2),"%"),"\r")
	OUT[i,(ncol(OUT)-7):ncol(OUT)]<-as.character(get.sp(BEST[i,1]))
	if(is.factor(OUT[1,10]))
		break()
	}

OUT[which(is.na(OUT[,13])),13]<-"unidentified"

OUT
}
