\name{simplify.network}
\alias{simplify.network}
\title{Network showing modules as nodes}
\description{This function modifies node coordinates to allow a clearer depiction of complex networks. Nodes are moved along the axis connecting the original position to the module centroid. The magnitude of such movement is defined by user.}
\usage{
simplify.network(node.names=NA,modules=NA,coordinates=NA,network=NA,
shift = 0.5,max.lwd.edge =2,min.lwd.edge =1,max.vertex.size=4,
min.vertex.size=2,label.size=1/2.5,bgcol="white",main="")
}
\arguments{
\item{node.names}{a vector containing the names of nodes}
\item{modules}{a vector containing the module assigned to each node}
\item{coordinates}{a two columns matrix containing the X and Y coordinates of each node in the original network}
\item{network}{a matrix describing the original network. Can be either a 0/1 matrix or a weighted matrix. Row names must contain node names.}
\item{shift}{a numeric defining the magnitude of node shift,limited between 0 (coinciding with the original location) and 1 (coinciding with the module centroid location).}
\item{max.lwd.edge}{if shift=1,a numeric defining the line width for the maximum number of connections between modules}
\item{min.lwd.edge}{if shift=1,a numeric defining the line width for the minimum number of connections between modules}
\item{max.vertex.size}{if shift=1,a numeric defining the size of the node representing the largest module}
\item{min.vertex.size}{if shift=1,a numeric defining the size of the node representing the smallest module}
\item{label.size}{a numeric defining the size of node labels,referred to its particular node size}
\item{bgcol}{a vector of strings representing the background colour for each node}
\item{main}{a string,the title for the plot (no title by default)}
}
\details{If 'shift=1',all nodes belonging to a module are represented as a single node depicted in the module centroid. In that case,node size is proportional to the number of element in this module and edge widths are proportional to the number of connections found between modules.}
\author{A. J. Muñoz-Pajares}
\examples{

inputMatrix<-matrix(c(1,1,1,1.2,2,1,0.8,1,3,2,1.2,1,4,2,2,2.2,
 5,3,1.8,2,6,3,2.2,2,7,3,1.7,2.1,8,3,2.2,2.2),ncol=4,byrow=TRUE)
colnames(inputMatrix)<-c("node","module","x","y")

network<-matrix(c(1,1,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,
  1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,
  0,1,1,1,1),ncol=8)
colnames(network)<-c(1:8)
row.names(network)<-c(1:8)

i1<-0
simplify.network(node.names=inputMatrix[,1],modules=inputMatrix[,2],
coordinates=inputMatrix[,3:4],network=network,shift = i1,
bgcol=c("red","red","blue","blue","green","green","green","green"),
main=paste("shift=",i1))

i1<-0.5
simplify.network(node.names=inputMatrix[,1],modules=inputMatrix[,2],
coordinates=inputMatrix[,3:4],network=network,shift = i1,
bgcol=c("red","red","blue","blue","green","green","green","green"),
main=paste("shift=",i1))

i1<-1.0
simplify.network(node.names=inputMatrix[,1],modules=inputMatrix[,2],
coordinates=inputMatrix[,3:4],network=network,shift = i1,
bgcol=c("red","red","blue","blue","green","green","green","green"),
main=paste("shift=",i1))

network<-as.matrix(as.dist(matrix(sample(c(1,0),10000,replace=TRUE),ncol=100)))
inputMatrix<-matrix(nrow=100,ncol=4)
inputMatrix[,1]<-1:100
inputMatrix[,2]<-c(rep(1,30),rep(2,20),rep(3,20),rep(4,20),rep(5,10))
inputMatrix[,3]<-c(
sample(seq(-40,0,0.01),30,rep=TRUE),
sample(seq(-40,0,0.01),20,rep=TRUE),
sample(seq(0,40,0.01),20,rep=TRUE),
sample(seq(0,40,0.01),20,rep=TRUE),
sample(seq(-20,20,0.01),10,rep=TRUE))
inputMatrix[,4]<-c(
sample(seq(0,40,0.01),30,rep=TRUE),
sample(seq(-40,0,0.01),20,rep=TRUE),
sample(seq(0,40,0.01),20,rep=TRUE),
sample(seq(-40,0,0.01),20,rep=TRUE),
sample(seq(-20,20,0.01),10,rep=TRUE))
cols<-c("red","green","yellow","blue","turquoise")

simplify.network(node.names=inputMatrix[,1],network=network,shift=0,
coordinates=inputMatrix[,3:4],modules=inputMatrix[,2],bgcol=cols[inputMatrix[,2]])

simplify.network(node.names=inputMatrix[,1],network=network,shift=1,
coordinates=inputMatrix[,3:4],modules=inputMatrix[,2],bgcol=cols[inputMatrix[,2]])

}

