using Distributed
using Plots
using ProgressMeter
using Distributions
using Printf
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
    t, vals = update(t, vals)
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

function nextDiv(rng)
  1.0 + randn(rng)/10.0
end

function divideCell(sc,rng,pSN)
 wildtypes = vcat(fill(true,sc.vals[1]),fill(false,sc.vals[2]))
 d1 = rand(Bool,length(wildtypes))
 d2 = [!x for x in d1]
 daughter1 = wildtypes[d1]
 daughter2 = wildtypes[d2]
 born1 = (vals = [sum(daughter1),sum([!x for x in daughter1])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng))
 born2 = (vals = [sum(daughter2),sum([!x for x in daughter2])], t0 = sc.tdiv, tdiv = sc.tdiv + nextDiv(rng))

 fate = rand()
 pSS = (1.0 - pSN)/2.0
 if fate <= pSS # S -> 2S
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
@everywhere rng, update = prepareSimulation(target = ss, seed = nothing)

time = @elapsed resarr = pmap(lifespan,fill([ss,0],24*5))

function plotRes(res)
  plot(res.age/365.0,res.mut+res.wt,title="mtDNA single cell",lw=2,legend = false,xlabel="Age (y)",ylabel="Copy Number",ylims = (0,1.25*ss));
  plot!(res.age/365.0,res.mut);
  plot!(res.age/365.0,res.wt);
  plot!(size = (1200,1200));
end

nrow = 3
ncol = 4
for i = 1:Int(ceil(length(resarr)/(nrow*ncol)))
  parr = plot([plotRes(res) for res in resarr[((nrow*ncol)*(i-1) + 1):((nrow*ncol)*i)]] ...,layout=(nrow,ncol),aspect_ratio = 0.1)
  fname = @sprintf("test%02d.png",i)
  savefig(parr,fname)
end

https://groups.google.com/d/msg/julia-users/0cV6v-FJD7c/eQcxNKWRAgAJ
func = interpolate((r.age,),r.wt, Gridded(Linear()));

ages = collect(0:100) # range(0,100,length=101)
reports = [interpolate((r.age,),r.wt, Gridded(Linear()))(ages)


# stemcells = [
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng)),
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng)),
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng)),
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng)),
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng)),
# (vals = [400,0], t0 = 0.0, tdiv = nextDiv(rng))
# ]

# times = [0.0]
# results = []

# while (times[end] < 500.0) & (length(stemcells) > 0)
  # global stemcells
  # global results
  # global times
  # push!(results,copy(stemcells))
  # dividingIndex = argmin([x.tdiv for x in stemcells])
  # dividingCell = stemcells[dividingIndex]
  # deleteat!(stemcells,dividingIndex)
  # push!(times,dividingCell.tdiv)
  # append!(stemcells,divideCell(dividingCell,rng,0.75))
# end
