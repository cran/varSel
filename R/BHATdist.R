#' @title Bhattacharyya distance among classes
#' @author Michele Dalponte and Hans Ole Oerka
#' @description Bhattacharyya distance.
#' @param g A column vector of the lables. length(g) is equal to nrow(X).
#' @param X A dataframe of the features. ncol(X) is equal to the total number of features, and nrow(X) is equal to the number of avaialble training samples. nrow(X) is equal to length(g)
#' @return A list containing a matrix of the class combinations and a vector of the Bhattacharyya distances of all the class combinations.
#' @export BHATdist
#' @references Dalponte, M., Oerka, H.O., Gobakken, T., Gianelle, D. & Naesset, E. (2013). Tree Species Classification in Boreal Forests With Hyperspectral Data. IEEE Transactions on Geoscience and Remote Sensing, 51, 2632-2645.
###############################################################################

BHATdist <- function(g,X){

  X<-as.matrix(X)

  nfeat <- ncol(X)
  nclass <- length(unique(g))

  mu <- by(X,g,colMeans)

  Cov <- by(X,g,stats::cov)

  ncomb <- t(utils::combn(unique(g),2))
  Bhat <- c()
  for(j in 1:nrow(ncomb)){
    mu.i <- mu[[ncomb[j,1]]]
    cov.i <- Cov[[ncomb[j,1]]]
    mu.j <- mu[[ncomb[j,2]]]
    cov.j <- Cov[[ncomb[j,2]]]
    if(nfeat==1){
      Bhat[j]<-(1/8)*t(mu.i-mu.j) %*% (solve((cov.i+cov.j)/2)) %*% (mu.i-mu.j) + 0.5*log((((cov.i+cov.j)/2))/(sqrt(((cov.i))*((cov.j)))),base=exp(1))
    }else{
      Bhat[j]<-(1/8)*t(mu.i-mu.j) %*% (solve((cov.i+cov.j)/2)) %*% (mu.i-mu.j) + 0.5*log(det(((cov.i+cov.j)/2))/(sqrt((det(cov.i))*(det(cov.j)))),base=exp(1))
    }
  }

  return(list(classComb=ncomb,bhatdist=Bhat))

}
