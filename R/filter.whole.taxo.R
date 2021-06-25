filter.whole.taxo<-function(whole.taxo){
OUT<-whole.taxo
OUT2<-cbind(OUT,"filtered_kingdom"=NA,"filtered_phylum"=NA, "filtered_class"=NA,"filtered_order"=NA,"filtered_family"=NA, "filtered_genus"=NA,"filtered_species"=NA)

for (i in 1:nrow(OUT2))
	{
	cat(paste(round(100*i/nrow(OUT2),2),"%"),"\r")

	this.line<-OUT[i,]
	whole.taxonomy<-this.line[7:13]

	OUT2[i,(ncol(OUT2)-6):(ncol(OUT2))]<-whole.taxonomy
	# Si hay celdas en blanco le pongo nombre.anterior.conocido:incertiae
	if(any(whole.taxonomy==""))
		{
		sonBlancos<-which(whole.taxonomy=="")
#		sonNada<-which(whole.taxonomy=="unidentified")
#		ultimo.conocido<-(min(sonBlancos,sonNada)-1)
		ultimo.conocido<-(min(sonBlancos)-1)
		whole.taxonomy[sonBlancos]<-paste(whole.taxonomy[ultimo.conocido],"incertiae",sep=":")
		OUT2[i,(ncol(OUT2)-6):(ncol(OUT2))]<-whole.taxonomy
		}

	if(this.line$pident>=97)
		level.retain<-21
	if(this.line$pident<97 & this.line$pident>=90)
		level.retain<-20
	if(this.line$pident<90 & this.line$pident>=85)
		level.retain<-19
	if(this.line$pident<85 & this.line$pident>=80)
		level.retain<-18
	if(this.line$pident<80 & this.line$pident>=75)
		level.retain<-17
	if(this.line$pident<75 & this.line$pident>=0)
		level.retain<-16

	if(level.retain<20)
		OUT2[i,level.retain:20]<-"low_pident"
	}	

OUT2
}
