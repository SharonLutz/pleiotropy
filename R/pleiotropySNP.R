
pleiotropySNP<-function(X,Y,Ydist,Z=NULL,covariates=FALSE,nPerm=5000){
    ###################
    # Input
    ###################
    # X is the a vector of the SNP of interest
    # Y is a matrix of the traits of interest where the number of rows equal the number of subjects and the number of columns equal the number of traits.
    # Ydist is a vector that specifies the distribution of the each trait. For example, for one normally distributed trait and the second binary trait, then Ydist<-c("gaussian","binomial"). Other choices for Ydist can be found by looking at the choice of family for the glm function (i.e. ?glm).
    # Z is a matrix of covariates where the number of rows equal the number of subjects and the number of columns equal the number of covariates.
    # covariates=FALSE then the models will not be adjusted for covariates, covariates=TRUE then the model will be adjusted for covariates
    # nPerm is the number of permutations. Default is 5,000.
    
    ###################
    #Output
    ###################
    # cutoffPvalue is the p-value from the cut-off based permutation approach
    # hausdorffPvalue is the p-value from the Hausdorff based permutation approach

  # TRUE if y is a matrix, FALSE if y is not a matrix
  isYmatrix = F
  if(!is.null(nrow(Y)) & !is.null(ncol(Y))){isYmatrix=TRUE}
  
  # TRUE if x is a matrix, FALSE if x is not a matrix
  isXmatrix = F
  if(!is.null(nrow(X)) & !is.null(ncol(X))){isXmatrix=TRUE}
  
  # If Z is a matrix, format appropriately
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
  
# CHECK that X is a vector
if(is.vector(X)&!isXmatrix){
    
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
         observedP[yIndex]<-summary(glm(Y[,yIndex]~X+Z),family=paste(Ydist[yIndex]))$coef[2,4]
        }else{
            observedP[yIndex]<-summary(glm(Y[,yIndex]~X),family=paste(Ydist[yIndex]))$coef[2,4]
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
   
   if(floor(jk/500)==ceiling(jk/500)){print(paste("permutation",jk, "out of ",nPerm))}
   
set.seed(jk)
Xperm<-sample(X,n) #permute X
permutedP<-rep(0,ncol(Y))#vector of p-values for observed

#calculate the p-value for each trait with permuted X
for(yIndex in 1:length(Ydist)){
    if(covariates==TRUE){
        permutedP[yIndex]<-summary(glm(Y[,yIndex]~Xperm+Z),family=paste(Ydist[yIndex]))$coef[2,4]
    }else{
        permutedP[yIndex]<-summary(glm(Y[,yIndex]~Xperm),family=paste(Ydist[yIndex]))$coef[2,4]
    }
}


###########################
#Cut-off based approach
############################

cutV<-(permutedP>observedP)
cutV0<-0
for(cutIndex in 1:length(observedP)){
    if(cutV[cutIndex]==TRUE){cutV0<-cutV0+1}
}
if(length(cutV)==cutV0  ){matPerm[jk,"PermC"]<-1}

############################
#Hausdorff based approach
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
    print("X (ie. SNP) is not a vector")}

list(cutoffPvalue=permP,hausdorffPvalue=hausP)}#end of function













