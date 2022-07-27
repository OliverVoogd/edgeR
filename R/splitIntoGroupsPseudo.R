splitIntoGroupsPseudo<-function(pseudo,group,pair) 
# Written by Davis McCarthy, February 2009, idea suggested by Mark Robinson
# A function to extract the data for specified two groups from a matrix of pseudocounts pair <- levels(as.factor(pair))
{
	y1<-pseudo[,group==pair[1]]; if (is.vector(y1)) { y1<-matrix(y1,ncol=1) }
	y2<-pseudo[,group==pair[2]]; if (is.vector(y2)) { y2<-matrix(y2,ncol=1) }
	y<-list(y1=y1,y2=y2)
	y
}