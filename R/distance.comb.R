distance.comb <-
function(matrices=NA,alphas=NA,method="Corrected",saveFile=TRUE)
{
#require(ape)

nmat<-length(matrices)
mats<-list()
length(mats)<-nmat
CORmats<-list()
length(CORmats)<-nmat

for (i in 1:nmat)
{
mats[[i]]<-eval(parse(text=matrices[i]))
CORmats[[i]]<-mats[[i]]/(max(mats[[i]]))
}

MSUM<-mats[[1]]*alphas[[1]]

for (i in 2:nmat)
{
if(method=="Uncorrected")
	{
	MATS<-mats[[i]]*alphas[[i]]
	MSUM<-MSUM+MATS
	OUTuncor<-MSUM
		if(saveFile==TRUE)
		write.table(OUTuncor,file="Weighted_Matrix_Uncorrected")
	}

if(method=="Corrected")
	{
	MATS<-CORmats[[i]]*alphas[[i]]
	MSUM<-MSUM+MATS
	OUTcor<-MSUM
		if(saveFile==TRUE)
		write.table(OUTcor,file="Weighted_Matrix_Corrected")
	}
}

if(method=="Uncorrected")
OUT<-OUTuncor

if(method=="Corrected")
OUT<-OUTcor

print(OUT)

}
