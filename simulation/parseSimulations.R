library(rjson)
mutf = "mutarr_519883308.json"
totf = "totarr_519883308.json"
muts = readChar(mutf, file.info(mutf)$size)
tots = readChar(totf, file.info(totf)$size)
muts = gsub("null",-99,muts)
tots = gsub("null",-99,tots)
mutarr = as.matrix(do.call(cbind,rjson::fromJSON(muts)))
totarr = as.matrix(do.call(cbind,rjson::fromJSON(tots)))
mutarr[mutarr==-99]=NA
totarr[mutarr==-99]=NA
time = seq(0,2,length.out=dim(mutarr)[1])

plot(time,mutarr[,1],type="n",ylim=c(0,max(mutarr,na.rm=TRUE)),ylab="Population size")
for(i in 1:dim(mutarr)[2]){
  points(time,mutarr[,i],type="l")
}