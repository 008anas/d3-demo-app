require(randomForest)

argv=commandArgs(T)

plus=paste(argv[1], "_plus.txt", sep="")
minus=paste(argv[1], "_minus.txt", sep="")

paste(plus)

scores_p=read.table(plus, stringsAsFactors=F, header=F)
scores_m=read.table(minus, stringsAsFactors=F, header=F)

colnames(scores_p)=c("Pos", "Score_-45", "Score_Energy", "Base Penalties", "TG", "Score_pribnow", "Score -35", "GC", "Strand")
rownames(scores_p)=scores_p[,1]

colnames(scores_m)=c("Pos", "Score_-45", "Score_Energy", "Base Penalties", "TG", "Score_pribnow", "Score -35", "GC", "Strand")
rownames(scores_m)=scores_m[,1]

scores_p[is.infinite(scores_p[,2]),2]=100
scores_m[is.infinite(scores_m[,2]),2]=100

equality=function(vector){
  a=which(vector[1:10]==vector[11])
  if(length(a)>1)
    a=a[1]
  a=26-a
  return(a)
}

X=dim(scores_p)[1]-26

min35=cbind(scores_p[c(X:dim(scores_p)[1],1:(X-1)),7], scores_p[c((X+1):dim(scores_p)[1],1:X),7], scores_p[c((X+2):dim(scores_p)[1],1:(X+1)),7], scores_p[c((X+3):dim(scores_p)[1],1:(X+2)),7], scores_p[c((X+4):dim(scores_p)[1],1:(X+3)),7], scores_p[c((X+5):dim(scores_p)[1],1:(X+4)),7], scores_p[c((X+6):dim(scores_p)[1],1:(X+5)),7],scores_p[c((X+7):dim(scores_p)[1],1:(X+6)),7],scores_p[c((X+8):dim(scores_p)[1],1:(X+7)),7],scores_p[c((X+9):dim(scores_p)[1],1:(X+8)),7])
minus35=apply(min35,1,min)
minus35=cbind(min35, minus35)
distminus35=apply(minus35, 1, equality)
scores_p[,10]=minus35[,11]
scores_p[,11]=distminus35


min35=cbind(scores_m[c(X:dim(scores_m)[1],1:(X-1)),7], scores_m[c((X+1):dim(scores_m)[1],1:X),7], scores_m[c((X+2):dim(scores_m)[1],1:(X+1)),7], scores_m[c((X+3):dim(scores_m)[1],1:(X+2)),7], scores_m[c((X+4):dim(scores_m)[1],1:(X+3)),7], scores_m[c((X+5):dim(scores_m)[1],1:(X+4)),7], scores_m[c((X+6):dim(scores_m)[1],1:(X+5)),7],scores_m[c((X+7):dim(scores_m)[1],1:(X+6)),7],scores_m[c((X+8):dim(scores_m)[1],1:(X+7)),7],scores_m[c((X+9):dim(scores_m)[1],1:(X+8)),7])
minus35=apply(min35,1,min)
minus35=cbind(min35, minus35)
distminus35=apply(minus35, 1, equality)
scores_m[,10]=minus35[,11]
scores_m[,11]=distminus35



colnames(scores_p)=c("Pos", "Score_-45", "Score_Energy", "Base_penalties", "TG", "Score_Pribnow","IGNORE", "GC", "Strand", "Score_-35", "Dist_-35")
colnames(scores_m)=colnames(scores_p)
data_rf=get(load(argv[2]))
result_p=predict(data_rf, type="prob", newdata=scores_p[,c(2:6,8,10)])
result_m=predict(data_rf, type="prob", newdata=scores_m[,c(2:6,8,10)])
result_p=cbind(scores_p, result_p[,2])
result_m=cbind(scores_m, result_m[,2])
  
  #superCoiling6h=read.delim(argv[3], header=F, stringsAsFactors=F, sep="\t")
  #rownames(superCoiling6h)=superCoiling6h[,1]
  #scores_p=cbind(scores_p, superCoiling6h[as.character(scores_p[,1]),2])
  #scores_m=cbind(scores_m, superCoiling6h[as.character(scores_m[,1]),2])
  #colnames(scores_p)=c("Pos", "Score_-45", "Score_Energy", "Base_penalties", "TG", "Score_Pribnow","IGNORE", "GC", "Strand", "Score_-35", "Dist_-35", "Supercoiling")
  #colnames(scores_m)=colnames(scores_p)
  #data_rf=get(load("RandomForestSupercoiling.RData"))
  #if(argv[2]=="S" | argv[2]=="s"){
  #  result_p=predict(data_rf, type="prob", newdata=scores_p[,c(2:6,8,10,12)])
  #  result_m=predict(data_rf, type="prob", newdata=scores_m[,c(2:6,8,10,12)])
  #  result_p=cbind(scores_p, result_p[,2])
  #  result_m=cbind(scores_m, result_m[,2])
  #} else {
  #  result_p=predict(data_rf, newdata=scores_p[,c(2:6,8,10,12)])
  #  result_m=predict(data_rf, newdata=scores_m[,c(2:6,8,10,12)])
  #  result_p=cbind(scores_p, result_p)
  #  result_m=cbind(scores_m, result_m)
  #}



result_p=result_p[,-7]
result_m=result_m[,-7]

write.table(result_p, file=paste(argv[1], "_final_plus.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
write.table(result_m, file=paste(argv[1], "_final_minus.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
