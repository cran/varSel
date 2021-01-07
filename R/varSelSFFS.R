#' @title Sequential Forward Floating Selection using Jeffries-Matusita Distance
#' @author Michele Dalponte and Hans Ole Oerka
#' @description Feature selection using the Sequential Forward Floating Selection search strategy and the Jeffries-Matusita distance.
#' @param g A column vector of the lables. length(g) is equal to nrow(X).
#' @param X A dataframe of the features. ncol(X) is equal to the total number of features, and nrow(X) is equal to the number of avaialble training samples. nrow(X) is equal to length(g)
#' @param strategy string indicating the multiclass strategy to adopt: 'minimum' or 'mean'.
#' @param n integer indicating the number of features to select. The algorithm will stop at n+1 features selected.
#' @return A list containing a vector of the JM distances on the individual bands, a matrix with the set of features selected, and a vector containing the distances for each feature set from 1 to N-1, where N is equal to ncol(X).
#' @export varSelSFFS
#' @references Dalponte, M., Oerka, H.O., Gobakken, T., Gianelle, D. & Naesset, E. (2013). Tree Species Classification in Boreal Forests With Hyperspectral Data. IEEE Transactions on Geoscience and Remote Sensing, 51, 2632-2645.
#' @examples
#' \dontrun{
#' data(dat)
#'
#' se<-varSelSFFS(g=dat$SP,X=dat[,c(1:65)],strategy="mean",n=4)
#' summary(se)
#'
#' }
###############################################################################


varSelSFFS <- function(g,X,strategy="mean",n=ncol(X)){

  if (strategy %in% c("minimum","mean")){

    if (n!=ncol(X)){
      STOP <- n+1
    }else{
      STOP = n
    }

    nfeat = ncol(X)

    nclass <- length(unique(g))

    mu <- by(X,g,colMeans)
    Cov <- by(X,g,stats::cov)

    #Jeffries-Matusita (JM) distance
    JM1<-function(fe){
      f<-fe
      ncomb <- t(utils::combn(1:nclass,2))
      Bhat <- c()
      jm <- c()
      for(j in 1:nrow(ncomb)){
        mu.i <- mu[[ncomb[j,1]]]
        cov.i <- Cov[[ncomb[j,1]]]
        mu.j <- mu[[ncomb[j,2]]]
        cov.j <- Cov[[ncomb[j,2]]]
        Bhat[j]<-(1/8)*t(mu.i[f]-mu.j[f]) %*% (solve((cov.i[f,f]+cov.j[f,f])/2)) %*% (mu.i[f]-mu.j[f]) + 0.5*log((((cov.i[f,f]+cov.j[f,f])/2))/(sqrt(((cov.i[f,f]))*((cov.j[f,f])))),base=exp(1))
        jm[j] <- sqrt(2*(1-exp(-Bhat[j])))
      }
      jm.dist<-c(min(jm),mean(jm))
      return(jm.dist)
    }

    JM<-function(fe){
      f<-fe
      ncomb <- t(utils::combn(1:nclass,2))
      Bhat <- c()
      jm <- c()
      for(j in 1:nrow(ncomb)){
        mu.i <- mu[[ncomb[j,1]]]
        cov.i <- Cov[[ncomb[j,1]]]
        mu.j <- mu[[ncomb[j,2]]]
        cov.j <- Cov[[ncomb[j,2]]]
        Bhat[j]<-(1/8)*t(mu.i[f]-mu.j[f]) %*% (solve((cov.i[f,f]+cov.j[f,f])/2)) %*% (mu.i[f]-mu.j[f]) + 0.5*log(det(((cov.i[f,f]+cov.j[f,f])/2))/(sqrt((det(cov.i[f,f]))*(det(cov.j[f,f])))),base=exp(1))
        jm[j] <- sqrt(2*(1-exp(-Bhat[j])))
      }
      jm.dist<-c(min(jm),mean(jm))
      return(jm.dist)
    }


    #Backword function
    backword<-function(fs){

      distance<-1:length(fs)

      for (j in 1:length(fs)){
        dist<-tryCatch(JM(fs[-j])[1],error=function(e){NA})
        if (!is.na(dist)){
          distance[j]<-dist
        }else{
          return<--1
        }
      }
      if (which.max(distance)==length(fs)){
        return<--1
      }else{
        return<-which.max(distance)
      }
    }


    features<-matrix(STOP,STOP,data=NA)
    distances<-matrix(1,STOP,data=NA)


    ###########################################################################################
    #STRART
    ###########################################################################################
    #computing JM for single features
    JMsingle<-matrix(nfeat,2,data=0)
    for(f in 1:nfeat){
      JMsingle[f,] <- JM1(f)
    }

    if(strategy == "minimum"){ref <- 1}

    if(strategy == "mean"){ref <- 2}

    features[1,1]<-which.max(JMsingle[,ref])
    distances[1] <- JMsingle[which.max(JMsingle[,ref]),1]

    print("Feature to select: 1")
    print(paste("JM distance: ",distances[1]))
    print("Features selected: ")
    print(features[1,1])

    if (STOP>1){

      #Two features
      i<-2
      C<-utils::combn(nfeat,2)
      dist<-apply(C,2,FUN=JM)

      FS<-C[,which.max(dist[ref,])]

      features[2,c(1,2)]<-C[,which.max(dist[ref,])]
      distances[2]<-dist[ref,which.max(dist[ref,])]

      print("Feature to select: 2")
      print(paste("JM distance: ",dist[ref,which.max(dist[ref,])]))
      print("Features selected: ")
      print(C[,which.max(dist[ref,])])

      if (STOP>2){

        #More than two features
        i<-2

        while (i<(nfeat-1) & (distances[i]<=round(sqrt(2),4)) & (i<STOP)){
          C<-matrix(i+1,(nfeat-i),data=NA)
          C[c(1:i),c(1:(nfeat-i))]<-replicate((nfeat-i),FS)
          ff<-1:nfeat
          ff<-ff[-FS]
          C[i+1,]<-ff

          dist<-tryCatch(apply(C,2,FUN=JM),error=function(e){NA})
          if (is.na(dist)[1]){
            print("ERROR: singularity of the covariance matrix")
            return(list(JMsingle=JMsingle[,ref],features=features,distances=distances,strategy=strategy, distance = "Jeffries-Matusita distance"))
          }

          FS<-C[,which.max(dist[ref,])]

          bw<-backword(FS)

          if (bw!=-1){

            FS<-FS[-bw]
            features[i,c(1:(i))]<-FS
            distTmp<-tryCatch(JM(FS),error=function(e){NA})
            if (is.na(distTmp)[1]){
              print("ERROR: singularity of the covariance matrix")
              return(list(JMsingle=JMsingle[,ref],features=features,distances=distances,strategy=strategy, distance = "Jeffries-Matusita distance"))
            }else{
              distances[i]<-distTmp[ref]
            }

            print(paste("Feature to select: ",i))
            print(paste("JM distance: ",distances[i]))
            print("Features selected: ")
            print(features[i,c(1:(i))])

          }else{

            features[i+1,c(1:(i+1))]<-C[,which.max(dist[ref,])]
            distances[i+1]<-dist[ref,which.max(dist[ref,])]
            i<-i+1
            print(paste("Feature to select: ",i))
            print(paste("JM distance: ",distances[i]))
            print("Features selected: ")
            print(features[i,c(1:(i))])

          }

        }

        return(list(JMsingle=JMsingle[,ref],features=features,distances=distances,strategy=strategy, distance = "Jeffries-Matusita distance"))

      }

    }

  }else{
    print("ERROR: wrong multiclass strategy")
  }
}

