##comun: A community table with samples as rows and taxa as columns. 
Com.niche <- function(comun, stats=TRUE){
    require(spaa)
    comun<-comun[,colSums(comun)>0]
    B<-niche.width(comun,method="levins")
    B_com<-1:nrow(comun)
    for(i in 1:nrow(comun)){
       a<-comun[i,]
       a<-a[a>0]
       B_com[i]<-mean(as.numeric(B[,names(a)]))
    }
    return(B_com)
}