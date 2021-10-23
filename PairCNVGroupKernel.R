KmerPairSim <- function(x,y,mersize){     ##x,y are matrices with 4 columns and kmersize rows
  
  Pairsim=list()
  for(m in 1:mersize){
    if (length(interval_intersection(interval(x[m,1],x[m,2]),interval(y[m,1],y[m,2])))==0) { Pairsim[m]=0 }
    if (length(interval_intersection(interval(x[m,1],x[m,2]),interval(y[m,1],y[m,2])))!=0){
      
      JacIndex= interval_measure(interval_intersection(interval(x[m,1],x[m,2]),interval(y[m,1],y[m,2])))/interval_measure(interval_union(interval(x[m,1],x[m,2]),interval(y[m,1],y[m,2])))   ##Contribution of length
      TypeCont= (((x[m,3]==y[[m,3]])+1)/2)     ##Contribution of type
      DosCont= 1/2^abs((abs(2-x[m,4])-abs(2-y[m,4])))     ##contribution of dosage
      Pairsim[m]= (JacIndex*TypeCont*DosCont)/mersize
      
    }
    
  }
  
  return((Reduce("+",Pairsim)))
}