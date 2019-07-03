using Turing
using StatsPlots
using Distributions
using Interpolations
using Random

function det_logistic(t, x0, r, K)
  (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))
end

function full_stoch_logistic(x0,r,K)
  eventno = K-x0
  unifs = Random.rand(eventno)
  clist = collect((x0+1):K)
  dts = -log.(1 .- unifs) ./ (r .* clist .* (1 .- clist ./ K))
  DataFrame(t=[0;cumsum(dts)], c=[x0;clist])
end

function stoch_logistic(t,x0,r,K)
  df = full_stoch_logistic(x0,r,K)
  interp = LinearInterpolation(df.t,df.c)
  interp(t)
end

x0 = 1
r = 3
K = 100
sigma = 5
logistic = stoch_logistic

# Generate some fake data
times =  collect(0:0.5:4)
data = [logistic(t,x0,r,K)+rand(Normal(0, sigma)) for t in times]
plot(times,data,seriestype=:scatter,title="Logistic model of population dynamics",lw=2,legend = false,xlabel="Time (d)",ylabel="Population size")

# Define a model with measurement error
@model fitlog(times,data) = begin
  s ~ InverseGamma(2,3)
  K ~ Uniform(0,500)
  r ~ Uniform(0,20)
  x0 ~ Uniform(0,10)
  for i in eachindex(times)
    data[i] ~ Normal(logistic(times[i],x0,r,K),s)
  end
end

chn = sample(fitlog(times,data), PG(1000,2000))
#chn = sample(fitlog(times,data), HMC(10000,0.1,5))
#chn = sample(fitlog(times,data), NUTS(1000,200,0.5))


# Summarise results (currently requires the master branch from MCMCChains)
describe(chn)

# Plot and save results
p = plot(chn)
savefig("gdemo-plot.png")

post = DataFrame(chn)
N = size(post,1)
prior = DataFrame(x0 = rand(Uniform(0,10),N), r = rand(Uniform(0,10),N), K = rand(Uniform(0,500),N), s = rand(InverseGamma(2,3),N))

chdat = Array{Float64}(undef,size(post,1),size(post,2),2)
chdat[:,:,1]=convert(Matrix,prior)
chdat[:,:,2]=Array(chn)
chn2 = Chains(chdat,chn.name_map.parameters)
plot(chn2)