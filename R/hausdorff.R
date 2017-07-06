hausdorff <-
function(Yob,Yperm){
    #compare the observed set of p-values to the permuted
    # using the hausdorff function
    # max(sup_x inf_y d(x,y), sup_y inf_x d(x,y))
    
    
    ###################
    #Input
    ###################
    #Yob is a vector of numbers (i.e. p-values for the observed traits with the SNP of interest)
    #Yperm is a vector of numbers (i.e. p-values for the permuted traits with the SNP of interest)
    
    ###################
    #Output
    ###################
    # hausdorff metric for the 2 sets
    
    minVop<-minVpo<-matrix(0,nrow=length(Yob),ncol=1)
    colnames(minVop)<-colnames(minVpo)<-c("Manh")
    
    for(jh1 in 1:length(Yob)){
        Ys1<-Yob[jh1]
        minV1<-matrix(0,nrow=length(Yperm),ncol=1)
        colnames(minV1)<-c("Manh")
        
        for(jh2 in 1:length(Yperm)){
            Ys2<-Yperm[jh2]
            minV1[jh2,"Manh"]<-abs(Ys1-Ys2)
        }
        minVop[jh1,"Manh"]<-min(minV1[,"Manh"])
    }
    
    for(jh2 in 1:length(Yperm)){
        Ys1<-Yperm[jh2]
        minV2<-matrix(0,nrow=length(Yob),ncol=1)
        colnames(minV2)<-c("Manh")
        
        for(jh1 in 1:length(Yob)){
            Ys2<-Yob[jh1]
            
            minV2[jh1,"Manh"]<-abs(Ys1-Ys2)
        }
        minVpo[jh2,"Manh"]<-min(minV2[,"Manh"])
    }
    
    Omanh<-max(c(max(minVop[,"Manh"]),max(minVpo[,"Manh"])))
    return(Omanh)}
