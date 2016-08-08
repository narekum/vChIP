args <- commandArgs(TRUE)
arg_length=length(args)

#fwdfile=read.table(file="high.fwd.zmatrix")
#revfile=read.table(file="high.rev.zmatrix")
#before=10
#after=20
fwdfile=args[1]
revfile=args[2]
bothfile=args[3]

before=as.numeric(args[4])
after=as.numeric(args[5])

positions=seq(-before,after)

readFile=function (file) {
	data=read.table(file=file)
	return(data)
}

getCutMatrix=function ( file) {
	noofcolumns=ncol(file)
	cutsmatrix=file[,seq(7,noofcolumns)]
	return (cutsmatrix)
}

getScores=function ( cutsmatrix ) {
	scores=as.vector(colMeans(cutsmatrix))
	return (scores)
}

fwddata=readFile(fwdfile)
revdata=readFile(revfile)
bothdata=readFile(bothfile)

fwdcutmatrix=getCutMatrix(fwddata)
revcutmatrix=getCutMatrix(revdata)
bothcutmatrix=getCutMatrix(bothdata)

fwdscores=getScores(fwdcutmatrix)
revscores=getScores(revcutmatrix)
bothscores=getScores(bothcutmatrix)

minvalue=min(fwdscores,revscores,bothscores)
maxvalue=max(fwdscores,revscores,bothscores)

pdf(file=paste(bothfile,".pdf",sep=""))
plot (positions,bothscores,type="l",ylim=c(minvalue,maxvalue))
lines(positions,fwdscores,col="blue")
lines(positions,revscores,col="red")
dev.off()

scores=cbind(positions,bothscores,fwdscores,revscores)
write.table(scores,file=paste(bothfile,".scores.dat",sep=""),quote=F,row.names = F,col.names = F)

