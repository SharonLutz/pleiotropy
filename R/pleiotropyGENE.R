pleiotropyGENE<-function(X,Y,Ydist,Z=NULL,covariates=FALSE,nPerm=1000){

library(SKAT)

###################
# Input
###################
# X is a matrix of rare and/or common variants in a region where the number of rows of the matrix equal the number of subjects and the number of columns equal the number of SNPs in the region.
# Y is a matrix of the traits of interest where the number of rows equal the number of subjects and the number of columns equal the number of traits.
# Ydist is a vector that specifies the distribution of the each trait. For example, for one normally distributed trait and the second binary trait, then Ydist<-c("C","D"). These choices are specified by the SKAT function.
# Z is a matrix of covariates where the number of rows equal the number of subjects and the number of columns equal the number of covariates.
 # covariates=FALSE then the models will not be adjusted for covariates, covariates=TRUE then the model will be adjusted for covariates.
# nPerm is the number of permutations. Default is 1,000.

###################
# Output
###################
# cutoffPvalue is the p-value from the cut-off based permutation approach
# hausdorffPvalue is the p-value from the Hausdorff based permutation approach

###################
# Warning
###################
# SKAT R package is needed to run this function. Make sure the SKAT package has been installed first.
  
  isXmatrix = F
  if(!is.null(nrow(X)) & !is.null(ncol(X))){
    # Format X to be a matrix
    temp <- NULL
    temp <- matrix(0,nrow=nrow(X),ncol=ncol(X))
    for(pp in 1:ncol(X)){temp[,pp]<-X[,pp]}
    isXmatrix=TRUE
    X <- temp
  }
  
  isYmatrix = F
  if(!is.null(nrow(Y)) & !is.null(ncol(Y))){
    # Format Y to be a matrix
    temp <- NULL
    temp <- matrix(0,nrow=nrow(Y),ncol=ncol(Y))
    for(pp in 1:ncol(Y)){temp[,pp]<-Y[,pp]}
    isYmatrix=TRUE
    Y <- temp
  }
  
  # Format Z to be a matrix (even if it is a vector) one covariate case
  if(!is.null(Z)){
    temp <- NULL
    if(!is.null(nrow(Z)) & !is.null(ncol(Z))){
      # Z is a matrix, format it to be recognized by SKAT
      temp <- matrix(0,nrow=nrow(Z),ncol=ncol(Z))
      for(pp in 1:ncol(Z)){temp[,pp]<-Z[,pp]}
    }else{
      # Z is a vector format to be a matrix
      temp <- matrix(0,nrow=length(Z),ncol=1)
      for(pp in 1:length(Z)){temp[pp,]<-Z[pp]}
    }
    # Set Z to be matrix
    Z <- temp
  }

  # CHECK that X is a matrix
if(isXmatrix){
    
    # CHECK that Y is a matrix
    if(isYmatrix){
        
        # CHECK that there are at least two traits
        if(ncol(Y)>1){
            
            # CHECK that Yval vector matches the dimnesions of the Y matrix
            if(ncol(Y)==length(Ydist)){
                
                observedP<-rep(0,ncol(Y))#vector of p-values for observed
                
                #calculate the observed p-value for each trait with X
                for(yIndex in 1:length(Ydist)){
                    if(covariates==TRUE){
                        obj1Z<-SKAT_Null_Model(Y[,yIndex] ~ Z, out_type=paste(Ydist[yIndex]))
                        observedP[yIndex]<-SKAT(X, obj1Z)$p.value
                    }else{
                         obj1nZ<-SKAT_Null_Model(Y[,yIndex] ~ 1, out_type=paste(Ydist[yIndex]))
                    observedP[yIndex]<-SKAT(X, obj1nZ)$p.value
                    }
                }

############################
# Permutation Loop
############################

#store permuted p-values for the 2 approaches
matPerm<-matrix(0,nrow=nPerm,ncol=2)
colnames(matPerm)<-c("PermC","HausP")

n<-nrow(Y)

#cycle through nPerm permutations
for(jk in 1:nPerm){
    
    if(floor(jk/100)==ceiling(jk/100)){print(paste("permutation",jk, "out of ",nPerm))}
    
    set.seed(jk)
    Xid<-sample(c(1:nrow(X)),n)
    
    Xperm<-X[Xid,] #permute X
    permutedP<-rep(0,ncol(Y))#vector of p-values for observed
    
    #calculate the p-value for each trait with permuted X
    for(yIndex in 1:length(Ydist)){
        if(covariates==TRUE){
            obj1Z<-SKAT_Null_Model(Y[,yIndex] ~ Z, out_type=paste(Ydist[yIndex]))
            permutedP[yIndex]<-SKAT(Xperm, obj1Z)$p.value
        }else{
            obj1nZ<-SKAT_Null_Model(Y[,yIndex] ~ 1, out_type=paste(Ydist[yIndex]))
            permutedP[yIndex]<-SKAT(Xperm, obj1nZ)$p.value
        }
    }


###########################
# Cut-off based approach
############################

cutV<-(permutedP>observedP)
cutV0<-0
for(cutIndex in 1:length(observedP)){
    if(cutV[cutIndex]==TRUE){cutV0<-cutV0+1}
}
if(length(cutV)==cutV0  ){matPerm[jk,"PermC"]<-1}

############################
# Hausdorff based approach
############################

originP<-rep(0,length(observedP))
OBmanh<-hausdorff(observedP,originP) #Hausdorff of observed to origin
Pmanh<-hausdorff(observedP,permutedP) #Hausdorff of observed to permuted
if(Pmanh<OBmanh){matPerm[jk,"HausP"]<-1}

############################
}#end of perm loop

permP<-1-sum(matPerm[,"PermC"])/nPerm #p-value using permutation cut-off method
hausP<-(sum(matPerm[,"HausP"])/nPerm)*(ncol(Y)^ncol(Y)+2) #p-value using hausdorff method
if(hausP>1){hausP<-1}

            }else{
                hausP<-NA
                permP<-NA
                print("The number of columns of Y does not match the length of Ydist")}
        }else{
            hausP<-NA
            permP<-NA
            print("There are not at least 2 traits (i.e. ncol(Y) are not greater than 1)")}
    }else{
        hausP<-NA
        permP<-NA
        print("Y is not a matrix")}
}else{
    hausP<-NA
    permP<-NA
    print("X is not a matrix")}

    list(cutoffPvalue=permP,hausdorffPvalue=hausP)}#end of function

