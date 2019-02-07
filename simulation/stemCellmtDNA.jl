using Distributed
using Plots
using ProgressMeter
using Distributions
using Statistics
using Printf
using Interpolations
import JSON
import Dates
@everywhere using Random

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

@everywhere function lifespan(vals)
  ages = collect(0.0:10.0:365.0*100.0)
  Nrecs = length(ages)
  res = (age = ages, wt = fill(0,Nrecs), mut = fill(0,Nrecs))
  t = 0.0

  #@showprogress 1 "Simulating..." for rec in 1:Nrecs
  for rec in 1:Nrecs
   t,vals = simulate(t,vals,res.age[rec])
   res.age[rec] = t
   res.wt[rec] = vals[1]
   res.mut[rec] = vals[2]  
  end
  return(res)
end

@everywhere function simulate(t, vals, tmax)
  while t < tmax
    if sum(vals) > 0
      t, vals = update(t, vals)
	else
	  t, vals = tmax, vals
	end
  end
  return(t,vals)
end

@everywhere function simulateSC(sc)
  tdiv, vals = simulate(sc.t0, sc.vals, sc.tdiv)
  return((vals = vals, t0 = sc.t0, tdiv = tdiv))
end

@everywhere function simulateSC(sc)
  return((vals = [2*x for x in sc.vals], t0 = sc.t0, tdiv = sc.tdiv))
end

function nextDiv(rng; mu=5, sigma=0.1)
  return(max(0,rand(rng,Normal(mu,sigma))))
end

function divideCell(sc, rng, pSN; N = 7, targ = 7, delta = 2, mu = 10, sigma = 1)
 # Randomly assign division type, except when overall stem cell population (N) has become too small (l.t. targ - delta) or too large (g.t. targ + delta)
 wildtypes = vcat(fill(true,sc.vals[1]),fill(false,sc.vals[2]))
 d1 = rand(Bool,length(wildtypes))
 d2 = [!x for x in d1]
 daughter1 = wildtypes[d1]
 daughter2 = wildtypes[d2]
 born1 = (vals = [sum(daughter1),sum([!x for x in daughter1])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng, mu, sigma))
 born2 = (vals = [sum(daughter2),sum([!x for x in daughter2])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng, mu, sigma))

 fate = rand()
 pSS = (1.0 - pSN)/2.0
 if (N < targ - delta) # S -> 2S
   cell1 = simulateSC(born1)
   cell2 = simulateSC(born2)
   return([cell1,cell2])
 elseif (N > targ + delta) # S -> 2N
   return([])
 elseif fate <= pSS  # S -> 2S
   cell1 = simulateSC(born1)
   cell2 = simulateSC(born2)
   return([cell1,cell2])
 elseif fate > pSS + pSN # S -> 2N
   return([])
 else
   cell1 = simulateSC(born1)
   return([cell1]) # S -> S + N
 end
end

function plotRes(res, ymin = 0, ymax=1000)
  plot(res.age/365.0,res.mut+res.wt,title="mtDNA single cell",lw=2,legend = false,xlabel="Age (y)",ylabel="Copy Number",ylims = (ymin,ymax));
  plot!(res.age/365.0,res.mut);
  plot!(res.age/365.0,res.wt);
  plot!(size = (1200,1200));
end

function removeDefective(fracarr, thresh = 1.0)
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

function getQuant(x,q)
  x = x[findall(.!isnan.(x))]
  if length(x) > 0
    res = quantile(x,q)
  else
    res = NaN
  end
  return(res)
end


@everywhere ss = 200
@everywhere rng, update = prepareSimulation(rfwt = 0.25, dwt = 0.01, rfmut = 0.25, dmut = 0.01, m = 3e-5, target = ss, seed = nothing)
#@everywhere rng, update = prepareSimulation(rfwt = 1.0, dwt = 0.14, rfmut = 1.0, dmut = 0.14, m = 3e-5, target = ss, seed = nothing)

println("Starting simulations!")
time = @elapsed resarr = pmap(lifespan,fill([ss/2,ss/2],24*600*9))
println(time)
print('\a')

# https://groups.google.com/d/msg/julia-users/0cV6v-FJD7c/eQcxNKWRAgAJ
# func = interpolate((r.age,),r.wt, Gridded(Linear()));

ages = collect(range(0,100*365,length=101))
mutarr = hcat([interpolate((r.age,),r.mut, Gridded(Linear()))(ages) for r in resarr]...);
totarr = hcat([interpolate((r.age,),r.mut+r.wt, Gridded(Linear()))(ages) for r in resarr]...);

anint = rand(0:999999999)
fname = @sprintf("mutarr_%09d.json",anint)
open(@sprintf("mutarr_%09d.json",anint),"w") do f
    JSON.print(f, mutarr);
end

open(@sprintf("totarr_%09d.json",anint),"w") do f
    JSON.print(f, totarr);
end

nrow = 3
ncol = 4
ymax = maximum([maximum([maximum(res.mut),maximum(res.wt),maximum(res.mut+res.wt)]) for res in resarr])
if ymax <= 100
 for i = 1:Int(ceil(length(resarr)/(nrow*ncol)))
   parr = plot([plotRes(res,0,ymax) for res in resarr[((nrow*ncol)*(i-1) + 1):((nrow*ncol)*i)]] ...,layout=(nrow,ncol))
   fname = @sprintf("test%02d.png",i)
   savefig(parr,fname)
 end
end

mutfrac = mutarr./totarr;
#fracarr = mutfrac
fracarr = removeDefective(mutfrac, 0.85)

low = [getQuant(fracarr[i,:],0.1) for i in 1:length(ages)];
mid = [getQuant(fracarr[i,:],0.5) for i in 1:length(ages)];
hig = [getQuant(fracarr[i,:],0.9) for i in 1:length(ages)];

mlab = @sprintf("Number of cells: %i",length(resarr))

quantplot = plot(ages/365.0,low,lw=1,legend=false,xlabel="Age (y)", ylabel="Mutation load", ylims = (0,1),title = mlab)
plot!(ages/365.0,mid,lw=2)
plot!(ages/365.0,hig,lw=1)
savefig(quantplot, "quantplot.png")

dynplot = plot(legend=false,xlabel="Age (y)", ylabel="Mutation load", ylims = (0,1),title = mlab)
plot!(ages/365.0, fracarr[:,1:1000], linealpha=0.1, linecolour=:black)
savefig(dynplot, "dynplot.png")

function overThresh(thresh, fracarr)
  [sum([x >= thresh for x in fracarr[i,:]])/length(fracarr[i,:]) for i in 1:size(fracarr)[1]]
end

threshplot = plot(legend=false,xlabel="Age (y)", ylabel="Fraction exceeding threshold", ylims = (0,1),title = mlab)
plot!(ages/365.0, overThresh(0.8,fracarr), linecolour=:black)
savefig(threshplot, "threshplot.png")
print('\a')

mu = 10;
sigma = 1;
stemcells = [
(vals = [100,100], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma)),
(vals = [100,100], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma)),
(vals = [100,010], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma)),
(vals = [100,100], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma)),
(vals = [100,100], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma)),
(vals = [100,100], t0 = 0.0, tdiv = nextDiv(rng, mu, sigma))
]

times = [0.0]
results = []

while (times[end] < 500.0) & (length(stemcells) > 0)
  global stemcells
  global results
  global times
  N = length(stemcells)
  push!(results,copy(stemcells))
  dividingIndex = argmin([x.tdiv for x in stemcells])
  dividingCell = stemcells[dividingIndex]
  deleteat!(stemcells,dividingIndex)
  push!(times,dividingCell.tdiv)
  append!(stemcells,divideCell(dividingCell, rng, 0.75; N = N, targ = 7, delta = 2, mu = mu, sigma = sigma))
end
