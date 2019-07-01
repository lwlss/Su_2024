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

dtsu = read.delim("../data/MouseData\ added number\ of mice\ 190219.csv",sep=",")
dmusc = read.delim("../data/CONNOR2_TG.csv",sep=",")

dcrypt = dtsu[dtsu$Tissue=="Epithelium",]
dcrypt$Mouse = as.integer(gsub("B","",dcrypt$Original.ID))

common = sort(intersect(unique(dcrypt$Mouse),unique(dmusc$Mouse)))

bw = 5.0
pch = "-"
datacol = rgb(1,0,0,0.35)
cex = 2
#bw = "nrd0"

musc_young = density(dmusc$MutationLoad[dmusc$Age==10],na.rm=TRUE, from = 0, to = 100)
musc_old = density(dmusc$MutationLoad[dmusc$Age==49],na.rm=TRUE, from = 0, to = 100)
simmusc_young = density(100*mutarr[which.min(abs(time/7 - 10)),]/totarr[which.min(abs(time/7 - 10)),],na.rm=TRUE, from = 0, to = 100,bw=bw)
simmusc_old = density(100*mutarr[which.min(abs(time/7 - 49)),]/totarr[which.min(abs(time/7 - 49)),],na.rm=TRUE, from = 0, to = 100,bw=bw)
muscmaxd = max(c(musc_young$y,musc_old$y,simmusc_young$y,simmusc_old$y))

crypt_young = density(dcrypt$MutationLoad[dcrypt$Age==10],na.rm=TRUE, from = 0, to = 100,bw=bw)
crypt_old = density(dcrypt$MutationLoad[dcrypt$Age==50],na.rm=TRUE, from = 0, to = 100,bw=bw)
simcrypt_young = density(100*mutarrcrypt[which.min(abs(time/7 - 10)),]/totarrcrypt[which.min(abs(time/7 - 10)),],na.rm=TRUE, from = 0, to = 100,bw=bw)
simcrypt_old = density(100*mutarrcrypt[which.min(abs(time/7 - 50)),]/totarrcrypt[which.min(abs(time/7 - 50)),],na.rm=TRUE, from = 0, to = 100,bw=bw)

cryptmaxd = max(c(crypt_young$y,crypt_old$y,simcrypt_young$y,simcrypt_old$y))


alpha = 0.2
#png("CompareDistributions.png", width = 2000, height=2000,pointsize=50)
pdf("CompareDistributions.pdf", useDingbats=FALSE)
op = par(mfrow=c(2,2),mar=c(4, 4, 3, 1) + 0.1)

plot(time/7,mutarr[,1],type="n",ylim=c(0,100),ylab="Mutation load (%)",main="Post-mitotic\nSkeletal muscle",xlab = "Age (weeks)")
for(i in 1:min(50,dim(mutarr)[2])){
  points(time/7,100*mutarr[,i]/totarr[,i],type="l",col=rgb(0,0,1,alpha ))
}
ages = dmusc$Age[dmusc$Age>=49]
points(dmusc$Age[dmusc$Age==10],dmusc$MutationLoad[dmusc$Age==10],pch=pch,col=datacol,cex=cex)
points(rep(mean(unique(ages)),length(ages)),dmusc$MutationLoad[dmusc$Age>=49],pch=pch,col=datacol,cex=cex)

plot(time/7,100*mutarr[,1],type="n",ylim=c(0,100),ylab="Mutation load (%)",main="Mitotic\nEpithelial crypt stem cells",xlab = "Age (weeks)")
for(i in 1:min(50,dim(mutarrcrypt)[2])){
  points(time/7,100*mutarrcrypt[,i]/totarrcrypt[,i],type="l",col=rgb(0,0,1,alpha ))
}
points(dcrypt$Age[dcrypt$Age==10],dcrypt$MutationLoad[dcrypt$Age==10],pch=pch,col=datacol,cex=cex)
points(dcrypt$Age[dcrypt$Age==50],dcrypt$MutationLoad[dcrypt$Age==50],pch=pch,col=datacol,cex=cex)

lwd = 2

plot(musc_young,xlim=c(0,100),ylim=c(0,muscmaxd),lwd=lwd,main="Post-mitotic\nSkeletal muscle",xlab="Mutation load (%)")
points(musc_old,type="l",col="red",lwd=lwd)
points(simmusc_young,type="l",lty=2,lwd=lwd)
points(simmusc_old,type="l",lty=2,col="red",lwd=lwd)

plot(crypt_young,xlim=c(0,100),ylim=c(0,cryptmaxd),lwd=lwd,main="Mitotic\nEpithelial crypt stem cells",xlab="Mutation load (%)")
points(crypt_old,type="l",col="red",lwd=lwd)
points(simcrypt_young,type="l",lty=2,lwd=lwd)
points(simcrypt_old,type="l",lty=2,col="red",lwd=lwd)

legend("topleft",c("young","old","experiment","simulation"),col=c("black","red","grey","grey"),lwd=lwd,lty=c(1,1,1,2),cex=0.75)

par(op)
dev.off()

png("CompareDistributions2.png", width = 2000, height=1000,pointsize=50)
op = par(mfrow=c(1,2),mar=c(4, 4, 3, 1) + 0.1)

plot(time/7,mutarr[,1],type="n",ylim=c(0,100),ylab="Mutation load (%)",main="Post-mitotic\nSkeletal muscle",xlab = "Age (weeks)")
for(i in 1:min(150,dim(mutarr)[2])){
  points(time/7,100*mutarr[,i]/totarr[,i],type="l",col=rgb(0,0,1,alpha ))
}
ages = dmusc$Age[dmusc$Age>=49]
points(dmusc$Age[dmusc$Age==10],dmusc$MutationLoad[dmusc$Age==10],pch=pch,col=datacol,cex=cex)
points(rep(mean(unique(ages)),length(ages)),dmusc$MutationLoad[dmusc$Age>=49],pch=pch,col=datacol,cex=cex)


plot(musc_young,xlim=c(0,100),ylim=c(0,muscmaxd),lwd=lwd,main="Post-mitotic\nSkeletal muscle",xlab="Mutation load (%)")
points(musc_old,type="l",col="red",lwd=lwd)
points(simmusc_young,type="l",lty=2,lwd=lwd)
points(simmusc_old,type="l",lty=2,col="red",lwd=lwd)

legend("topleft",c("young","old","experiment","simulation"),col=c("black","red","grey","grey"),lwd=lwd,lty=c(1,1,1,2),cex=0.75)

par(op)
dev.off()


pdf("Individuals.pdf", useDingbats=FALSE)

# Split figures
for (mouse in common){

op = par(mfrow=c(1,2),mar=c(4, 4, 3, 1) + 0.1,oma=c(0,0,0,0))

plot(time/7,mutarr[,1],type="n",ylim=c(0,100),ylab="Mutation load (%)",main="Post-mitotic\nSkeletal muscle",xlab = "Age (weeks)",cex.lab=1.5,cex.axis=1.5,cex.main=1.25)
for(i in 1:min(50,dim(mutarr)[2])){
  points(time/7,100*mutarr[,i]/totarr[,i],type="l",col=rgb(0,0,1,alpha ))
}
ages = dmusc$Age[(dmusc$Age>=49)&(dmusc$Mouse==mouse)]
points(dmusc$Age[(dmusc$Age==10)&(dmusc$Mouse==mouse)],dmusc$MutationLoad[(dmusc$Age==10)&(dmusc$Mouse==mouse)],pch=pch,col=datacol,cex=cex)
points(rep(mean(unique(ages)),length(ages)),dmusc$MutationLoad[(dmusc$Age>=49)&(dmusc$Mouse==mouse)],pch=pch,col=datacol,cex=cex)

plot(time/7,100*mutarr[,1],type="n",ylim=c(0,100),ylab="Mutation load (%)",main="Mitotic\nEpithelial crypt stem cells",xlab = "Age (weeks)",cex.lab=1.5,cex.axis=1.5,cex.main=1.25)
for(i in 1:min(50,dim(mutarrcrypt)[2])){
  points(time/7,100*mutarrcrypt[,i]/totarrcrypt[,i],type="l",col=rgb(0,0,1,alpha ))
}
points(dcrypt$Age[(dcrypt$Age==10)&(dcrypt$Mouse==mouse)],dcrypt$MutationLoad[(dcrypt$Age==10)&(dcrypt$Mouse==mouse)],pch=pch,col=datacol,cex=cex)
points(dcrypt$Age[(dcrypt$Age==50)&(dcrypt$Mouse==mouse)],dcrypt$MutationLoad[(dcrypt$Age==50)&(dcrypt$Mouse==mouse)],pch=pch,col=datacol,cex=cex)

par(op)
title(paste("Mouse ID:",mouse),outer=TRUE,line=-1)

}

dev.off()






# Counts for figure legend
length(dmusc$MutationLoad[dmusc$Age==10])
length(dmusc$MutationLoad[dmusc$Age>=49])
length(dcrypt$MutationLoad[dcrypt$Age==10])
length(dcrypt$MutationLoad[dcrypt$Age==50])

stime = time/7

# Calculations for Results section
round(mean(dmusc$MutationLoad[dmusc$Age==10],na.rm=TRUE))
sd(dmusc$MutationLoad[dmusc$Age==10],na.rm=TRUE)

round(mean(dmusc$MutationLoad[dmusc$Age>=49],na.rm=TRUE))
sd(dmusc$MutationLoad[dmusc$Age>=49],na.rm=TRUE)

stime[10]
round(mean(100*mutarr[10,]/totarr[10,],na.rm=TRUE))
sd(100*mutarr[10,]/totarr[10,],na.rm=TRUE)

stime[48]
round(mean(100*mutarr[48,]/totarr[48,],na.rm=TRUE))
sd(100*mutarr[48,]/totarr[48,],na.rm=TRUE)

round(mean(dcrypt$MutationLoad[dcrypt$Age==10],na.rm=TRUE))
sd(dcrypt$MutationLoad[dcrypt$Age==10],na.rm=TRUE)

round(mean(dcrypt$MutationLoad[dcrypt$Age==50],na.rm=TRUE))
sd(dcrypt$MutationLoad[dcrypt$Age==50],na.rm=TRUE)

stime[10]
round(mean(100*mutarrcrypt[10,]/totarrcrypt[10,],na.rm=TRUE))
sd(100*mutarrcrypt[10,]/totarrcrypt[10,],na.rm=TRUE)

stime[49]
round(mean(100*mutarrcrypt[49,]/totarrcrypt[49,],na.rm=TRUE))
sd(100*mutarrcrypt[49,]/totarrcrypt[49,],na.rm=TRUE)
