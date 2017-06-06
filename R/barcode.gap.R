#setwd("/home/ajesusmp/Dropbox/Barcode_COI/2015_10_thresholds/04_2015_10_20")
#inter<-read.table(file="out.interNew.txt")
#intra<-read.table(file="out.intraNew.txt")


barcode.gap<-function(summary=NULL,stat.intra="max",stat.inter="min",
xlab=NULL, ylab=NULL, legend=TRUE, lab.nodes="nogap")
{
if(length(summary)!=0)
{
INTRA<-as.data.frame(summary$Intraspecific)
INTER<-as.data.frame(summary$Interspecific)

if(stat.intra=="max")
Intra<-INTRA$Max.
if(stat.intra=="min")
Intra<-INTRA$Min.
if(stat.intra=="median")
Intra<-INTRA$Median.
if(stat.intra=="mean")
Intra<-INTRA$Mean

if(stat.inter=="max")
Inter<-INTER$Max.
if(stat.inter=="min")
Inter<-INTER$Min.
if(stat.inter=="median")
Inter<-INTER$Median.
if(stat.inter=="mean")
Inter<-INTER$Mean

inter.names<-intra.names<-row.names(summary$Intraspecific)
if (length(xlab)==0)
xlab=paste(stat.intra,"intraspecific distances")
if (length(ylab)==0)
ylab=paste(stat.inter,"interspecific distances")
}

		T1<-unique(sort(as.matrix(Intra)))[2]
		T2<-T1+0.1*T1
		T3<-unique(sort(Inter))[2]
		T4<-T3+0.1*T3
		T5<-max(sort(Intra))
		T6<-max(sort(Inter))
		T7<-T5*0.8
		T8<-T6*0.8
		T9<-T5*1.05

	max.intra.plot<-max(Intra)
	max.inter.plot<-max(Inter)

if(legend==TRUE)
	{
	dev.new(width=10.5,height=6.993110)
	layout(matrix(rep(c(1,1,1,2),4),ncol=4,byrow=TRUE))
	par(mar=c(5.5,5.5,1.8,0.5))
		plot(Intra,Inter,pch=18,cex.axis=1.8,cex.lab=1.8, xlab=xlab, ylab=ylab,type="n",xlim=c(0,T9))
	
		polygon(x=c(-1,max(T5,T6),2*max(T5,T6)),y=c(-1,max(T5,T6),0),col="lightgrey")
		points(Intra,Inter,pch=18,cex=2.5)
		abline(a=0,b=1)
		box(lty=1)

		outs.nogap<-cbind(Intra[which(Intra>Inter)],Inter[which(Intra>Inter)])
		row.names(outs.nogap)<- intra.names[which(Intra>Inter)]
	
		if(length(outs.nogap)!=0)
			{
			if(is.null(row.names(outs.nogap))) row.names(outs.nogap)<-paste("sp",c(1:nrow(outs.nogap)))
			colnames(outs.nogap)<-c("Intra","Inter")
	
			aux1.nogap<-sapply(row.names(outs.nogap),function(x){strsplit(x," ")})
			nombres.nogap<-paste(sapply(aux1.nogap,function(x){substr(x[1],1,1)}),". ",sapply(aux1.nogap,function(x){x[2]}),sep="")
			}

		outs.gap<-cbind(Intra[which(Intra<Inter)],Inter[which(Intra<Inter)])
		row.names(outs.gap)<- intra.names[which(Intra<Inter)]

		if(length(outs.gap)!=0)
			{
			if(is.null(row.names(outs.gap))) row.names(outs.gap)<-paste("sp",c(1:nrow(outs.gap)))
			colnames(outs.gap)<-c("Intra","Inter")

			aux1.gap<-sapply(row.names(outs.gap),function(x){strsplit(x," ")})
			nombres.gap<-paste(sapply(aux1.gap,function(x){substr(x[1],1,1)}),". ",sapply(aux1.gap,function(x){x[2]}),sep="")
			}

		if(lab.nodes=="nogap"|lab.nodes=="all")
		if(length(outs.nogap)!=0)
		text(outs.nogap,labels=c(1:nrow(outs.nogap)),pos=4,cex=1.5)

		if(lab.nodes=="gap"|lab.nodes=="all")
		if(length(outs.gap)!=0)
		text(outs.gap,labels=c(1:nrow(outs.gap)),pos=4,cex=1.5)

	par(mar=c(1.8,0,1.5,0.2))
	plot(c(1,1),type="n",axes=FALSE,xlab="",ylab="")

	# añadir otro plot finito al lado izq y poner la leyenda de colores y labels.
	legend(x="top",pch=22,col="black",pt.bg=c("white","lightgray"),legend=c("Barcode Gap","No Barcode Gap"),cex=1.6)

	if(lab.nodes=="nogap"|lab.nodes=="all")
		ifelse(length(outs.nogap)!=0,{
		legend(x="left",pch="",legend=paste(1:nrow(outs.nogap),row.names(outs.nogap)),cex=1.5,xjust=0.2,bty="n")},{
		legend(x="left",pch="",legend="All species show\n barcode gap",cex=1.5,xjust=0.2,bty="n")})

	if(lab.nodes=="gap"|lab.nodes=="all")
		ifelse(length(outs.gap)!=0,{
		legend(x="left",pch="",legend=paste(1:nrow(outs.gap),row.names(outs.gap)),cex=1.5,xjust=0.2,bty="n")},{
		legend(x="left",pch="",legend="All species lack\n barcode gap",cex=1.5,xjust=0.2,bty="n")})
	}

if(legend==FALSE)
	{
	plot(Intra,Inter,pch=18,xlab=xlab, ylab=ylab,type="n",xlim=c(0,T9))
	
	polygon(x=c(-1,max(T5,T6),2*max(T5,T6)),y=c(-1,max(T5,T6),0),col="lightgrey")
	points(Intra,Inter,pch=18)
	abline(a=0,b=1)
	box(lty=1)

	outs.nogap<-cbind(Intra[which(Intra>Inter)],Inter[which(Intra>Inter)])
	row.names(outs.nogap)<- intra.names[which(Intra>Inter)]
	
		if(length(outs.nogap)!=0)
			{
			if(is.null(row.names(outs.nogap))) row.names(outs.nogap)<-paste("sp",c(1:nrow(outs.nogap)))
			colnames(outs.nogap)<-c("Intra","Inter")
	
			aux1.nogap<-sapply(row.names(outs.nogap),function(x){strsplit(x," ")})
			nombres.nogap<-paste(sapply(aux1.nogap,function(x){substr(x[1],1,1)}),". ",sapply(aux1.nogap,function(x){x[2]}),sep="")
			}

		outs.gap<-cbind(Intra[which(Intra<Inter)],Inter[which(Intra<Inter)])
		row.names(outs.gap)<- intra.names[which(Intra<Inter)]

		if(length(outs.gap)!=0)
			{
			if(is.null(row.names(outs.gap))) row.names(outs.gap)<-paste("sp",c(1:nrow(outs.gap)))
			colnames(outs.gap)<-c("Intra","Inter")

			aux1.gap<-sapply(row.names(outs.gap),function(x){strsplit(x," ")})
			nombres.gap<-paste(sapply(aux1.gap,function(x){substr(x[1],1,1)}),". ",sapply(aux1.gap,function(x){x[2]}),sep="")
			}

	}

OUT<-list()
OUT[[1]]<-outs.nogap
OUT[[2]]<-outs.gap
names(OUT)<-c("no.barcode.gap","barcode.gap")
OUT
}

#barcode.gap(inter,intra)
#barcode.gap(inter$Min.,intra$Max.,names=row.names(intra))
#barcode.gap(as.matrix(inter$Min.),as.matrix(intra$Max.),names=row.names(intra))

####






#plot(rep(1,2*length(tabla[which(tabla$Ratio>10),1])),c(tabla[which(tabla$Ratio>10),1],tabla[which(tabla$Ratio>10),2]),ylim=c(0,0.4))

# The gap is between the max of intra and the min of inter:
#  gap10<-c(max(tabla[which(tabla$Ratio>10),1]),min(tabla[which(tabla$Ratio>10),2]))
#  gap4<-c(max(tabla[which(tabla$Ratio>4),1]),min(tabla[which(tabla$Ratio>4),2]))

# I estimate the mean point in this gap:

#threshold10<-diff(gap10)/2+min(gap10)
#threshold4<-diff(gap4)/2+min(gap4)


