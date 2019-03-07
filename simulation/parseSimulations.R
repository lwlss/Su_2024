library(rjson)
parseJSON = function(fname){
  bigstring = readChar(fname, file.info(fname)$size)
  bigstring = gsub("null",-99,bigstring)
  bigarr = as.matrix(do.call(cbind,rjson::fromJSON(bigstring)))
  bigarr[bigarr == -99] = NA
  return(bigarr)
}

flist = list.files()
flist = flist[grepl("totarr_crypt_",flist)]
flist = flist[order(file.info(flist)$mtime,decreasing=TRUE)]
rint = substring(flist[1],14,22)

mutarr = parseJSON(sprintf("mutarr_%s.json",rint))
totarr = parseJSON(sprintf("totarr_%s.json",rint))
mutarrcrypt = parseJSON(sprintf("mutarr_crypt_%s.json",rint))
totarrcrypt = parseJSON(sprintf("totarr_crypt_%s.json",rint))


time = seq(5,2*365,length.out=dim(mutarr)[1])

dcrypt = read.delim("../data/MouseData\ added number\ of mice\ 190219.csv",sep=",")
dmusc = read.delim("../data/CONNOR1_TG.csv",sep=",")

bw = 0.05
#bw = "nrd0"

musc_young = density(dmusc$MutationLoad[dmusc$Age==10]/100,na.rm=TRUE, from = 0, to = 1)
musc_old = density(dmusc$MutationLoad[dmusc$Age==49]/100,na.rm=TRUE, from = 0, to = 1)
simmusc_young = density(mutarr[which.min(abs(time/7 - 10)),]/totarr[which.min(abs(time/7 - 10)),],na.rm=TRUE, from = 0, to = 1,bw=bw)
simmusc_old = density(mutarr[which.min(abs(time/7 - 49)),]/totarr[which.min(abs(time/7 - 49)),],na.rm=TRUE, from = 0, to = 1,bw=bw)
muscmaxd = max(c(musc_young$y,musc_old$y,simmusc_young$y,simmusc_old$y))

crypt_young = density(dcrypt$MutationLoad[dcrypt$Age==10]/100,na.rm=TRUE, from = 0, to = 1,bw=bw)
crypt_old = density(dcrypt$MutationLoad[dcrypt$Age==50]/100,na.rm=TRUE, from = 0, to = 1,bw=bw)
simcrypt_young = density(mutarrcrypt[which.min(abs(time/7 - 10)),]/totarrcrypt[which.min(abs(time/7 - 10)),],na.rm=TRUE, from = 0, to = 1,bw=bw)
simcrypt_old = density(mutarrcrypt[which.min(abs(time/7 - 50)),]/totarrcrypt[which.min(abs(time/7 - 50)),],na.rm=TRUE, from = 0, to = 1,bw=bw)

cryptmaxd = max(c(crypt_young$y,crypt_old$y,simcrypt_young$y,simcrypt_old$y))

png("CompareDistributions.png", width = 3000, height=2000,pointsize=50)
op = par(mfrow=c(2,2))

plot(time/7,mutarr[,1],type="n",ylim=c(0,1),ylab="Mutation load",main="Post-mitotic\nSkeletal muscle",xlab = "Age (d)")
for(i in 1:min(50,dim(mutarr)[2])){
  points(time/7,mutarr[,i]/totarr[,i],type="l",col=rgb(0,0,1,0.3))
}
points(dmusc$Age[dmusc$Age==10],dmusc$MutationLoad[dmusc$Age==10]/100,pch=16,col=rgb(0,0,0,0.2))
points(dmusc$Age[dmusc$Age==49],dmusc$MutationLoad[dmusc$Age==49]/100,pch=16,col=rgb(1,0,0,0.2))

plot(time/7,mutarr[,1],type="n",ylim=c(0,1),ylab="Mutation load",main="Mitotic\nEpithelial crypt stem cells",xlab = "Age (d)")
for(i in 1:min(50,dim(mutarrcrypt)[2])){
  points(time/7,mutarrcrypt[,i]/totarrcrypt[,i],type="l",col=rgb(0,0,1,0.3))
}
points(dcrypt$Age[dcrypt$Age==10],dcrypt$MutationLoad[dcrypt$Age==10]/100,pch=16,col=rgb(0,0,0,0.2))
points(dcrypt$Age[dcrypt$Age==50],dcrypt$MutationLoad[dcrypt$Age==50]/100,pch=16,col=rgb(1,0,0,0.2))

lwd = 2

plot(musc_young,xlim=c(0,1),ylim=c(0,muscmaxd),lwd=lwd,main="Post-mitotic\nSkeletal muscle",xlab="Mutation load")
points(musc_old,type="l",col="red",lwd=lwd)
points(simmusc_young,type="l",lty=2,lwd=lwd)
points(simmusc_old,type="l",lty=2,col="red",lwd=lwd)

plot(crypt_young,xlim=c(0,1),ylim=c(0,cryptmaxd),lwd=lwd,main="Mitotic\nEpithelial crypt stem cells",xlab="Mutation load")
points(crypt_old,type="l",col="red",lwd=lwd)
points(simcrypt_young,type="l",lty=2,lwd=lwd)
points(simcrypt_old,type="l",lty=2,col="red",lwd=lwd)

legend("topleft",c("young","old","simulation","experiment"),col=c("black","red","grey","grey"),lwd=lwd,lty=c(1,1,1,2))

par(op)
dev.off()

png("CompareDistributions2.png", width = 3000, height=2000,pointsize=50)
op = par(mfrow=c(1,2))

plot(time/7,mutarr[,1],type="n",ylim=c(0,1),ylab="Mutation load",main="Post-mitotic\nSkeletal muscle",xlab = "Age (d)")
for(i in 1:min(150,dim(mutarr)[2])){
  points(time/7,mutarr[,i]/totarr[,i],type="l",col=rgb(0,0,1,0.3))
}
points(dmusc$Age[dmusc$Age==10],dmusc$MutationLoad[dmusc$Age==10]/100,pch=16,col=rgb(0,0,0,0.2))
points(dmusc$Age[dmusc$Age==49],dmusc$MutationLoad[dmusc$Age==49]/100,pch=16,col=rgb(1,0,0,0.2))

lwd = 3

plot(musc_young,xlim=c(0,1),ylim=c(0,muscmaxd),lwd=lwd,main="Post-mitotic\nSkeletal muscle",xlab="Mutation load")
points(musc_old,type="l",col="red",lwd=lwd)
points(simmusc_young,type="l",lty=2,lwd=lwd)
points(simmusc_old,type="l",lty=2,col="red",lwd=lwd)

legend("topleft",c("young","old","simulation","experiment"),col=c("black","red","grey","grey"),lwd=lwd,lty=c(1,1,1,2))

par(op)
dev.off()

