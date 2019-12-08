@everywhere using Distributed
using Plots
using Plots.PlotMeasures
using ProgressMeter
@everywhere using Distributions
@everywhere using Statistics
@everywhere using Printf
@everywhere using Interpolations
import JSON
import Dates
@everywhere using Random
using DataFrames
using CSV

@everywhere function prepareSimulation(;rfwt = 0.25, dwt = 0.14, rfmut = 0.25, dmut = 0.14, m = 3e-5, target = 3000, seed = nothing)
 if seed == nothing 
   seed = rand(0:999999999);
 end   
 rng = Random.MersenneTwister(seed);
 delta_wt = [1,-1,0,0,0];
 delta_mut = [0,0,1,1,-1];
 rwt = dwt * rfwt;
 rmut = dmut * rfmut;
 
 function update(t, vals)
   wt = vals[1];
   mut = vals[2];
   tot = wt + mut;
   Fwt = wt/tot;
   Fmut = mut/tot;
   haz = [((dwt - rwt)*target*Fwt + rwt*wt)*(1-m), dwt*wt, ((dwt - rwt)*target*Fwt + rwt*wt)*m, (dmut - rmut)*target*Fmut + rmut*mut, dmut*mut]; # control
   #haz = [(rwt*wt)*(1-m), dwt*wt, (rwt*wt)*m, rmut*mut, dmut*mut];
   hazfrac = cumsum(haz/sum(haz));
   index = minimum(findall(hazfrac .- Random.rand(rng) .> 0));
   t = t + (1.0/sum(haz))*(Random.randexp(rng));
   vals[1] = wt + delta_wt[index];
   vals[2] = mut + delta_mut[index];
   return((t,vals));
 end
 
 return(rng,update);
end

@everywhere function lifespan(vals, agemax; lethal_thresh = 1.0)
  ages = collect(0.0:1:(365.0*agemax))
  Nrecs = length(ages)
  res = (age = ages, wt = fill(0.0,Nrecs), mut = fill(0.0,Nrecs))
  t = 0.0

  #@showprogress 1 "Simulating..." for rec in 1:Nrecs
  for rec in 1:Nrecs
   t,vals = simulate(t,vals,res.age[rec], lethal_thresh = lethal_thresh)
   res.age[rec] = t
   res.wt[rec] = vals[1]
   res.mut[rec] = vals[2]  
  end
  return(res)
end

@everywhere function simulate(t, vals, tmax; lethal_thresh = 1.0)
  if isnan(vals[1]) || isnan(vals[2])
    return(tmax,vals)
  end
  while t < tmax
    N = sum(vals)
    if vals[2]/N <= lethal_thresh
      if sum(vals) > 0
        t, vals = update(t, vals)
	  else
	    t, vals = tmax, vals
	  end
	else 
	  t, vals = tmax, [NaN,NaN]
	end
  end
  return(t,vals)
end

function plotRes(res, ymin = 0, ymax=1000)
  plot(res.age/365.0,res.mut+res.wt,title="mtDNA single cell",lw=2,legend = false,xlabel="Age (y)",ylabel="Copy Number",ylims = (ymin,ymax));
  plot!(res.age/365.0,res.mut);
  plot!(res.age/365.0,res.wt);
  plot!(size = (1200,1200));
end

@everywhere function removeDefective(fracarr, thresh = 1.0)
  fracadj = copy(fracarr)
  for i in 1:size(fracarr)[2]
    ff = findfirst(fracarr[:,i].>thresh)
    if ff != nothing
      for j in ff:size(fracarr)[1]
       fracadj[j,i] = NaN
      end
    end
  end
  return(fracadj)
end  

@everywhere function getQuant(x,q)
  x = x[findall(.!isnan.(x))]
  if length(x) > 0
    if(q=="mean")
	  res = mean(x)
	else
	  res = quantile(x,q)
	end
  else
    res = NaN
  end
  return(res)
end

function initPop(source,N,ss)
  muts = [Int(round(x/100.0)) for x in sample(rng,source,N)*ss]
  wts = [ss-mut for mut in muts] 
  p0 = [[wt,mut] for (wt,mut) in zip(wts,muts)]
  return(p0)
end

dat = CSV.read("../data/MouseData.csv")
dat.Days = dat.Age * 7
dat2 = CSV.read("../data/CONNOR1_TG.csv")
dat2.Days = dat2.Age * 7
epi = dat[dat.Tissue .== "Epithelium",:]
mus = dat2[[x!==missing for x in dat2.MutationLoad],:]
inits = mus.MutationLoad[mus.Age.==10]

inits = 100*rand(Normal(0.9, 0.025), 1000)

@everywhere ss = 200
@everywhere agemax = 12 # years
@everywhere lethal_thresh = 0.975 #0.9
@everywhere mtDNA_deg = 0.0075 # 0.1
@everywhere rng, update = prepareSimulation(rfwt = 0.25, dwt = mtDNA_deg, rfmut = 0.25, dmut = mtDNA_deg, m = 3e-5, target = ss, seed = nothing)

println("Starting simulations!")
time = @elapsed resarr = pmap(x -> lifespan(x, agemax, lethal_thresh = lethal_thresh),initPop(inits,24*500,ss))
println(time)
print('\a')

# https://groups.google.com/d/msg/julia-users/0cV6v-FJD7c/eQcxNKWRAgAJ
# func = interpolate((r.age,),r.wt, Gridded(Linear()));

ages = collect(range(5,agemax*365,length=101))
mutarr = hcat([interpolate((r.age,),r.mut, Gridded(Linear()))(ages) for r in resarr]...);
totarr = hcat([interpolate((r.age,),r.mut+r.wt, Gridded(Linear()))(ages) for r in resarr]...);

anint = rand(0:999999999)
open(@sprintf("human_muscle_mutarr_%09d.json",anint),"w") do f
    JSON.print(f, mutarr);
end

open(@sprintf("human_muscle_totarr_%09d.json",anint),"w") do f
    JSON.print(f, totarr);
end

nrow = 3
ncol = 4
ymax = maximum([maximum([maximum(res.mut),maximum(res.wt),maximum(res.mut+res.wt)]) for res in resarr])
for i = 1:min(4,Int(ceil(length(resarr)/(nrow*ncol))))
   parr = plot([plotRes(res,0,ymax) for res in resarr[((nrow*ncol)*(i-1) + 1):((nrow*ncol)*i)]] ...,layout=(nrow,ncol))
   fname = @sprintf("human_muscle_mtDNAdyn%02d.png",i)
   savefig(parr,fname)
end


mutfrac = mutarr./totarr;
#fracarr = mutfrac
fracarr = removeDefective(mutfrac, 1.0)

low = [getQuant(fracarr[i,:],0.1) for i in 1:length(ages)];
mid = [getQuant(fracarr[i,:],0.5) for i in 1:length(ages)];
hig = [getQuant(fracarr[i,:],0.9) for i in 1:length(ages)];
mean_frac = [getQuant(fracarr[i,:],"mean") for i in 1:length(ages)];

getQuant(epi.MutationLoad[epi.Age.==10],0.1)
getQuant(epi.MutationLoad[epi.Age.==10],0.9)

mlab = @sprintf("Number of cells: %i",length(resarr))

cols = [:red,:green,:blue,:orange]
quantplot = plot(ages/365.0,low,lw=1,legend=false,xlabel="Age (y)", ylabel="Mutation load", ylims = (0,1),title = mlab, linecolor = :black, linestyle = :dash, bottom_margin = 10px)
plot!(ages/365.0,mid,lw=2, linecolor = :black)
plot!(ages/365.0,hig,lw=1, linecolor = :black, linestyle = :dash)
plot!(ages/365.0,mean_frac,lw=2, linecolor = :red)
#plot!(epi.Days/365.0,epi.MutationLoad/100,seriestype=:scatter,markercolor = [cols[m] for m in epi.Mouse],markerstrokecolor = false)
savefig(quantplot, "human_muscle_quantplot.png")

dynplot = plot(legend=false,xlabel="Age (y)", ylabel="Mutation load", ylims = (0,1),title = mlab)
plot!(ages/365.0, fracarr[:,1:min(size(fracarr)[2],1000)], linealpha=0.1, linecolour=:black)
savefig(dynplot, "human_muscle_dynplot.png")

function overThresh(thresh, fracarr)
  [sum([x >= thresh for x in fracarr[i,:]])/length(fracarr[i,:]) for i in 1:size(fracarr)[1]]
end

mkdir("threshplots")
defthreshes = [0.0:0.05:1.0;];
fsize = 18;
threshplot = plot(legend=false,xlabel="Time since first biopsy (y)", ylabel="Fraction exceeding threshold", ylims = (0,1),title = mlab,
xtickfontsize=fsize,ytickfontsize=fsize,xguidefontsize=fsize,yguidefontsize=fsize,legendfontsize=fsize,titlefontsize=fsize,size=(1000,750))
for thresh in defthreshes
  plot!(ages/365.0, overThresh(thresh,fracarr))
end
savefig(threshplot, "human_muscle_threshplot.png")