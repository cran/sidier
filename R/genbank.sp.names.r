genbank.sp.names<-function(sequences)
{
fun1<-function(x){paste(strsplit(x,split=" ")[[1]][c(2,3)],collapse=" ")}
new.names<-sapply(row.names(sequences),fun1)
names(new.names)<-c("")
row.names(sequences)<-new.names
sequences
}

