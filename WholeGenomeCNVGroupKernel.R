KmerTotal <-
  function(x,y){   ## x,y are list of matrices
    p=length(x)
    q=length(y)
    if(p*q==0){TotalSim=0}
    if(p==q&q>0){
      Pairsim=rep(0,p)
      for ( i in 1:p)
      {
        Pairsim[i]=KmerPairSim(x[[i]],y[[i]],1)      
      }
      TotalSim=sum(Pairsim)
    }
    if(p>q&q>0){
      SPairsim=rep(0,p-q+1)
      for ( l in 0:(p-q)){
        Pairsim=rep(0,q)
        for ( k in 1:q)
        {
          Pairsim[k]=KmerPairSim(x[[k+l]],y[[k]],1)     
          
        }
        SPairsim[(l+1)]=sum(Pairsim)
      }
      TotalSim=max(SPairsim)
    }
    
    if(q>p&p>0){
      SPairsim=rep(0,q-p+1)
      for ( l in 0:(q-p)){
        Pairsim=rep(0,p)
        for ( k in 1:p)
        {
          Pairsim[k]=KmerPairSim(x[[k]],y[[k+l]],1)       
          
        }
        SPairsim[(l+1)]=sum(Pairsim)
      }
      TotalSim=max(SPairsim)
    }
    return(TotalSim)
  }