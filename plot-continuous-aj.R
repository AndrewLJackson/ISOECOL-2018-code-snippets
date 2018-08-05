# my attempt to recover proportions from ilr output from mixsiar

# needs objects: jags.1, source and mix to be available

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Transform regression lines from ILR-space to p-space
e <- with( source, {
  
  # set up the e matrix
  e <- matrix(rep(0, n.sources * (n.sources - 1)),
              nrow = n.sources,
              ncol = (n.sources-1))
  
  for(i in 1:(n.sources-1)){
    
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i), 
                   -sqrt(i/(i+1)), 
                   rep(0,n.sources-i-1)))
    
    # divide by its sum
    e[,i] <- e[,i]/sum(e[,i])
  }
  
  return(e)
}) 


# the inverse here is the problem... its really really clunky
inverseIrl <- function(ilr_x, e) {
  
  cross.med <- array(data=NA,dim=c(nrow(ilr_x), n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p.med <- array(data=NA,dim=c(nrow(ilr_x), n.sources))              # dummy variable for inverse ILR calculation
  p.median <- array(data=NA,dim=c(nrow(ilr_x), n.sources))
  
  for(i in 1:nrow(ilr_x)){
    for(j in 1:(n.sources-1)){
      cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
    }
    for(src in 1:n.sources){
      tmp.p.med[i,src] <- prod(cross.med[i,src,]);
    }
    for(src in 1:n.sources){
      p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
    }
  }
  return(p.median)
}
  




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# get the columns ilr.contX and ilr.global out

idx.cont <- grep("ilr.cont", colnames(jags.1$BUGSoutput$sims.matrix)  )
ilr.cont <- jags.1$BUGSoutput$sims.matrix[, idx.cont]

idx.global <- grep("ilr.global", colnames(jags.1$BUGSoutput$sims.matrix) )
ilr.global <- jags.1$BUGSoutput$sims.matrix[, idx.global]
                   

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# define a range of value we want to evaluate our continuous variable over
xx <- seq(0, 7.5)

# calculate the ilr scaled y data
b0 <- ilr.global[1,]
b1 <- ilr.cont[1,]
yy <-  xx %*% t(b1) + rep(1,length(xx)) %*% t(b0)

pp <- inverseIrl(yy, e)

