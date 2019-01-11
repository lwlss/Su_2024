%% Random number generator from a user defined discrete probability 
% distribution

% x - vector of outcomes
% p - vector of outcome probabilities 
% ns - how many random numbers you need

function S = DiscSampVec2(x,p,ns) 
 
[~,idx] = histc(rand(1,ns),[0,cumsum(p)]); 
S = x(idx);
