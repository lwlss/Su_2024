using Distributed
using Plots
using ProgressMeter
using Distributions
@everywhere using Random

@everywhere function prepareSimulation(;rfwt = 0.25, dwt = 0.14, rfmut = 0.25, dmut = 0.14, m = 3e-5, target = 3000, seed = nothing)
 if seed == nothing 
   seed = rand(0:999999999);
 end   
 rng = Random.MersenneTwister(seed);
 delta_wt = [1,-1,-1,0,0];
 delta_mut = [0,0,1,1,-1];
 rwt = dwt * rfwt;
 rmut = dmut * rfmut;
 
 function update(t, vals)
   wt = vals[1];
   mut = vals[2];
   tot = wt + mut;
   Fwt = wt/tot;
   Fmut = mut/tot;
   haz = [((dwt - rwt)*target*Fwt + rwt*wt)*(1-m), dwt*wt, ((dwt - rwt)*target*Fwt + rwt*wt)*m, (dmut - rmut)*target*Fmut + rmut*mut, dmut*mut];
   hazfrac = cumsum(haz/sum(haz));
   index = minimum(findall(hazfrac .- Random.rand(rng) .> 0));
   t = t + (1.0/sum(haz))*(Random.randexp(rng));
   vals[1] = wt + delta_wt[index];
   vals[2] = mut + delta_mut[index];
   return((t,vals));
 end
 
 return(update);
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
    t, vals = update(t, vals)
  end
  return(t,vals)
end

@everywhere function simulateSC(sc)
  tdiv, vals = simulate(sc.t0, sc.vals, sc.tdiv)
  return((vals = vals, t0 = sc.t0, tdiv = tdiv))
end

function nextDiv()
  1.0 + randn()/10.0
end

function divideCell(sc,pSN = 0.75)
 wildtypes = vcat(fill(true,sc.vals[1]),fill(false,sc.vals[2]))
 d1 = rand(Bool,length(wildtypes))
 d2 = [!x for x in d1]
 daughter1 = wildtypes[d1]
 daughter2 = wildtypes[d2]
 born1 = (vals = [sum(daughter1),sum([!x for x in daughter1])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv())
 born2 = (vals = [sum(daughter2),sum([!x for x in daughter2])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv())

 fate = rand()
 pSS = (1.0 - pSN)/2.0
 if fate < pSS # S -> 2S
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

@everywhere ss = 400
@everywhere update = prepareSimulation(target = ss)

#time = @elapsed resarr = pmap(lifespan,fill([ss,0],4))

#res = resarr[1]
#plot(res.age/365.0,hcat(res.mut+res.wt,res.mut,res.wt),title="Stem cell mtDNA population dynamics",lw=2,legend = false,xlabel="Age (y)",ylabel="Copy Number")

stemcells = [
(vals = [400,0], t0 = 0.0, tdiv = nextDiv()),
(vals = [400,0], t0 = 0.0, tdiv = nextDiv()),
(vals = [400,0], t0 = 0.0, tdiv = nextDiv()),
(vals = [400,0], t0 = 0.0, tdiv = nextDiv()),
(vals = [400,0], t0 = 0.0, tdiv = nextDiv()),
(vals = [400,0], t0 = 0.0, tdiv = nextDiv())
]

times = [0.0]
results = []

while times[end] < 5.0
  append!(results,[stemcells])
  dividingIndex = argmin([x.tdiv for x in stemcells])
  dividingCell = stemcells[dividingIndex]
  deleteat!(stemcells,dividingIndex)
  append!(times,dividingCell.tdiv)
  append!(stemcells,divideCell(dividingCell))
end
  
  
  
  
  
  



