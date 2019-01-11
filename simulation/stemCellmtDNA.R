prepareSimulation = function(rfwt = 0.25, dwt = 0.14, rfmut = 0.25, dmut = 0.14, m = 3e-5, target = 3000){
 delta_wt = c(1,-1,-1,0,0)
 delta_mut = c(0,0,1,1,-1)
 rwt = dwt * rfwt
 rmut =dmut *rfmut

 update = function(vals){
   t = vals[1]
   wt = vals[2]
   mut = vals[3]
   tot = wt + mut
   Fwt = wt/tot
   Fmut = mut/tot
   haz = c(((dwt - rwt)*target*Fwt + rwt*wt)*(1-m), dwt*wt, ((dwt - rwt)*target*Fwt + rwt*wt)*m, (dmut - rmut)*target*Fmut + rmut*mut, dmut*mut)
   hazfrac = cumsum(haz/sum(haz))
   index = min(which(hazfrac-runif(1)>0))
   vals[1] = t + rexp(1,sum(haz))
   vals[2] = wt + delta_wt[index]
   vals[3] = mut + delta_mut[index]
   return(vals)
  }
}

#Most animal cells in culture appear to have
#approximately 1000–5000 molecules of the
#mitochondrial genome per cell
#https://doi.org/10.1016/S0168-9525(01)02238-7

ss = 3000

# Fig. 2 from Morales review
rdelta = 0.1
res = data.frame(age=seq(0,35,rdelta))
vals = c(0,100,0) # wt, mut

# Human lifespan
rdelta = 10
res = data.frame(age=seq(0,365*100,rdelta))
vals = c(0,3000,0) # wt, mut

res$wt = 0
res$mut = 0

Nrecs = length(res$age)
pb = txtProgressBar(min = 1, max = Nrecs, style = 3)
update = prepareSimulation(target = ss)
for(rec in 1:Nrecs){
 while(vals[1]<res$age[rec]) vals = update(vals)
 res[rec,] = vals
 setTxtProgressBar(pb, rec)
}
op = par(mfrow = c(1,2))
plot(res$age/365,res$wt+res$mut,type="l",xlab="Age (y)",ylab="# molecules",ylim=range(c(0,ss,res$wt,res$mut,res$wt+res$mut),na.rm=TRUE),lwd=2)
abline(h=ss,col="black",lty=2,lwd=2)
points(res$age/365,res$mut,type="l",col="red",lwd=2)
points(res$age/365,res$wt,type="l",col="blue",lwd=2)
#legend("left",c("wt","mutant","total"),col=c("blue","red","black"),lwd=2,bg="white")

plot(res$age/365,100*res$mut/(res$mut+res$wt),type="l",xlab="Age (y)",ylab="Mutation load (%)",ylim=c(0,100),lwd=2,col="red")
par(op)

# Need to think about how well controlled total mtDNA copy number is
# Basic Gillespie simulation has little to say about biologically plausible
# levels of variance in total copy number

f = function(x) ss/x
A = 10
b = (log(A)-1)/ss
f2 = function(x) A * exp(-b*x)
#curve(f,from=0,to=5000)
#curve(f2,from=0,to=5000)
