@everywhere using Distributed
using Plots
using ProgressMeter
@everywhere using Distributions
@everywhere using Statistics
@everywhere using Printf
@everywhere using Interpolations
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

@everywhere function lifespan(vals, agemax)
  ages = collect(0.0:1:(365.0*agemax))
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

#@everywhere function simulateSC(sc)
#  return((vals = [2*x for x in sc.vals], t0 = sc.t0, tdiv = sc.tdiv))
#end

@everywhere function nextDiv(rng; mu=5, sigma=0.1)
  return(max(0,rand(rng,Normal(mu,sigma))))
end

@everywhere function divideCell(sc, rng, pSN; N = 7, targ = 7, delta = 2, mu = 10, sigma = 1, m = 1e-5, before = true, defect_thresh = 1.0)
 # Randomly assign division type, except when overall stem cell population (N) has become too small (l.t. targ - delta) or too large (g.t. targ + delta)
 wildtypes = vcat(fill(true,sc.vals[1]),fill(false,sc.vals[2]))
 if before # Simulate duplication from N -> 2N before cell division, else, N -> N/2 and recovers towards target
  wt_copy = copy(wildtypes)
  for w in 1:length(wt_copy)
   if(wt_copy[w])
     if rand() < m 
	  wt_copy[w] = false
	 end
   end
  end
  wildtypes = vcat(wildtypes,wt_copy)
 end
 
 d1 = rand(Bool,length(wildtypes))
 d2 = [!x for x in d1]
 daughter1 = wildtypes[d1]
 daughter2 = wildtypes[d2]
 born1 = (vals = [sum(daughter1),sum([!x for x in daughter1])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng; mu=mu,sigma=sigma))
 born2 = (vals = [sum(daughter2),sum([!x for x in daughter2])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng; mu=mu,sigma=sigma))

 fate = rand()
 pSS = (1.0 - pSN)/2.0
 divprod =[]
 if (N < targ - delta) # S -> 2S
   cell1 = simulateSC(born1)
   cell2 = simulateSC(born2)
   divprod = [cell1,cell2]
 elseif (N > targ + delta) # S -> 2N
   divprod = []
 elseif fate <= pSS  # S -> 2S
   cell1 = simulateSC(born1)
   cell2 = simulateSC(born2)
   divprod = [cell1,cell2]
 elseif fate > pSS + pSN # S -> 2N
   divprod = []
 else
   cell1 = simulateSC(born1)
   divprod = [cell1] # S -> S + N
 end
 divprod = [c for c in divprod if c[2]/(c[1]+c[2]) <= thresh] # Simulate loss of daughters with biochemical deficiency
 return(divprod)
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
    res = quantile(x,q)
  else
    res = NaN
  end
  return(res)
end

@everywhere function simCrypt(;Nsc = 7, v0 = [100,100], t_days = 100, mu = 2, sigma = 0.1, m = 3e-5, defect_thresh = 1.0)
  stemcells = [(vals = v0, t0 = 0.0, tdiv = nextDiv(rng; mu=mu,sigma=sigma)) for i in 1:Nsc]

  times = [0.0]
  results = []

  while (times[end] < t_days) & (length(stemcells) > 0)
    N = length(stemcells)
    push!(results,copy(stemcells))
    dividingIndex = argmin([x.tdiv for x in stemcells])
    dividingCell = stemcells[dividingIndex]
    deleteat!(stemcells,dividingIndex)
    push!(times,dividingCell.tdiv)
    append!(stemcells,divideCell(dividingCell, rng, 0.75; N = N, targ = Nsc, delta = 2, mu = mu, sigma = sigma, m = m, before = true, defect_thresh = 0.85))
  end
  return(times,results)
end

@everywhere ss = 200
@everywhere agemax = 2 # years
@everywhere rng, update = prepareSimulation(rfwt = 0.25, dwt = 0.1, rfmut = 0.25, dmut = 0.1, m = 3e-5, target = ss, seed = nothing)
#@everywhere rng, update = prepareSimulation(rfwt = 1.0, dwt = 0.14, rfmut = 1.0, dmut = 0.14, m = 3e-5, target = ss, seed = nothing)

println("Starting simulations!")
time = @elapsed resarr = pmap(x -> lifespan(x, agemax),fill([60,140],24*50))
println(time)
print('\a')

# https://groups.google.com/d/msg/julia-users/0cV6v-FJD7c/eQcxNKWRAgAJ
# func = interpolate((r.age,),r.wt, Gridded(Linear()));

ages = collect(range(0,agemax*365,length=101))
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
plot!(ages/365.0, fracarr[:,1:min(size(fracarr)[2],1000)], linealpha=0.1, linecolour=:black)
savefig(dynplot, "dynplot.png")

function overThresh(thresh, fracarr)
  [sum([x >= thresh for x in fracarr[i,:]])/length(fracarr[i,:]) for i in 1:size(fracarr)[1]]
end

threshplot = plot(legend=false,xlabel="Age (y)", ylabel="Fraction exceeding threshold", ylims = (0,1),title = mlab)
plot!(ages/365.0, overThresh(0.8,fracarr), linecolour=:black)
savefig(threshplot, "threshplot.png")

timesc = @elapsed scarr = pmap(x -> simCrypt(Nsc = 7, v0 = [60,140], t_days = 365*agemax, mu = 2, sigma = 0.1, m = 3e-5, defect_thresh = 0.85),1:(24*12))
println(timesc)
print("\a")

ages = collect(range(5,agemax*365,length=101))
muts = []
tots = []
for i in 1:length(scarr)
  times, results = scarr[i]
  Ndyn = [(minimum([c.tdiv for c in sc]),length(sc)) for sc in results]
  mtDNA = [(minimum([c.tdiv for c in sc]),sum([c.vals[1] for c in sc]),sum([c.vals[2] for c in sc])) for sc in results]
  summ = [(t = minimum([c.tdiv for c in sc]), N = length(sc), wt = sum([c.vals[1] for c in sc]), mut = sum([c.vals[2] for c in sc])) for sc in results]
  mutint = interpolate(([s.t for s in summ],),[float(s.mut) for s in summ], Gridded(Linear()))(ages)
  totint = interpolate(([s.t for s in summ],),[float(s.wt + s.mut) for s in summ], Gridded(Linear()))(ages)
  append!(muts,[mutint])
  append!(tots,[totint])
  
  scparr = plot([
   plot([d[1]/365 for d in Ndyn],[d[2] for d in Ndyn], xlabel="Age (y)", ylabel="Number of stem cells", ylims=(0,maximum([d[2] for d in Ndyn])),legend=false),
   plot([d[1]/365 for d in mtDNA],[100*d[3]/(d[2]+d[3]) for d in mtDNA], xlabel="Age (y)", ylabel="Mutation load (%)", ylims=(0,100),legend=false)]
   ..., layout = (1,2) 
  )
  savefig(scparr,@sprintf("SC%02d.png",i))
end

fracs = hcat(muts...)./hcat(tots...)
#fracs = removeDefective(fracs,0.85)
fplot = plot(legend=false,xlabel="Age (y)", ylabel="Mutation load", ylims = (0,100))
plot!(ages/365,100*fracs,linealpha=0.5,linecolour=:black)

low = [getQuant(fracs[i,:],0.1) for i in 1:length(ages)];
mid = [getQuant(fracs[i,:],0.5) for i in 1:length(ages)];
hig = [getQuant(fracs[i,:],0.9) for i in 1:length(ages)];

mlab = @sprintf("Number of crypts: %i",size(fracs)[2])

quantplot = plot(ages/365.0,low,lw=1,legend=false,xlabel="Age (y)", ylabel="Mutation load (%)", ylims = (0,1),title = mlab)
plot!(ages/365.0,mid,lw=2)
plot!(ages/365.0,hig,lw=1)
#savefig(quantplot, "quantplot.png")


