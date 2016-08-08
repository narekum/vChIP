args<-commandArgs(TRUE)

stateFile=args[1]
predictionFile=args[2]



state=read.table(file=stateFile,sep="\t")
state=state[,1]

prediction=read.table(file=predictionFile,sep="\t")
prediction=prediction[,1]

library(pROC)

pdf(file=paste(stateFile,"roc.pdf",sep="."))
rocobj <- plot.roc(state,prediction,percent=TRUE,print.auc=TRUE )
dev.off()
