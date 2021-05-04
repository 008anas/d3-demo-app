argv=commandArgs(T)

promotersPlus=read.table(argv[1], header=T, stringsAsFactors=F, sep="\t")
promotersMinus=read.table(argv[2], header=T, stringsAsFactors=F, sep="\t")

predictedP=promotersPlus[promotersPlus[,dim(promotersPlus)[2]]>0.6,]
predictedM=promotersMinus[promotersMinus[,dim(promotersMinus)[2]]>0.6,]

predictedP=predictedP[order(predictedP[,1]),]

i=2
while(i<=dim(predictedP)[1]){
  if(abs(predictedP[i,1]-predictedP[(i-1),1])<18){
    a=which(predictedP[(i-1):i,11]==max(predictedP[(i-1):i,11]))
    if(length(a)>1){
      a=which(predictedP[(i-1):i,dim(predictedP)[2]]==max(predictedP[(i-1):i, dim(predictedP)[2]]))
    }
    if(a==1){
      predictedP[i,1]=NA
    } else{
      predictedP[(i-1),1]=NA
    }
    predictedP=na.omit(predictedP)
  } else{
    i=i+1
  }
}

predictedM=predictedM[order(predictedM[,1]),]

i=2
while(i<=dim(predictedM)[1]){
  if(abs(predictedM[i,1]-predictedM[(i-1),1])<18){
    a=which(predictedM[(i-1):i,11]==max(predictedM[(i-1):i,11]))
    if(length(a)>1){
      a=which(predictedM[(i-1):i,dim(predictedM)[2]]==max(predictedM[(i-1):i, dim(predictedM)[2]]))
    }
    if(a==1){
      predictedM[i,1]=NA
    } else{
      predictedM[(i-1),1]=NA
    }
    predictedM=na.omit(predictedM)
  } else{
    i=i+1
  }
}


colnames(predictedP)=colnames(predictedM)
predicted=rbind(predictedP, predictedM)
colnames(predicted[dim(predicted)[2]])="Score"



write.table(predicted, file=argv[3], sep="\t", col.names=T, quote=F, row.names=F )


# 
# predicted2=promotersPlus[promotersPlus[,11]>=0.75,]
# predictedM2=promotersMinus[promotersMinus[,11]>=0.75,]
# a=c()
# b=c()
# for(i in 25:(length(promotersPlus[,1])-25)){
#   a[i]=49+i
#   b[i]=sum(promotersPlus[i-24,11]:promotersPlus[i+25,11])
# }
# 
# calculate=function(i, vector){
#   a=prod(vector[i:(i+29)])
#   return(a)
# }
# 
# vac=unlist(lapply(1:816315,FUN=calculate, vector=promotersPlus[,11]))
# a=(1:816315)+15
