%% Random number generator from a user defined discrete probability 
% distribution

% x - vector of outcomes
% p - vector of outcome probabilities 
% ns - how many random numbers you need


function S = DiscSampVec3(x,p,ns)

global gg

S = zeros(gg.numDiv,gg.mtDNA*gg.initS);

for ii = 1 : gg.numDiv
    
    pVec = [1-p(ii),p(ii)];
    
    [~,idx] = histc(rand(1,ns),[0,cumsum(pVec)]);
    S(ii,:) = x(idx);
    
end

end
