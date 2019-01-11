function [MutatedSCAgeFission,CryptFisTime2] = mtDNACrypt_ConstantV2FCN_CF_Input(s3)

% mtDNACrypt Function for Constant and Increasing Mutation Rate
% - Crypt Fission - Single Crypts - FINAL

% The script is an amalgamation of the previous crypt model. It identifies
% that there are a certain number of mtDNA molecules residing within each
% stem cell of the crypt. With the evolution of stem cell divisions, the
% number of mutated mtDNA molecules evolves stochastically according to
% pre-determined probabilities. Also, with each additional mutated mtDNA
% molecule, the model determines which kind of mutation has developed
% according to probability data previously aquired. Therefore, this model
% is a more accurate representation of the processes that take place within
% the crypt and at the tissue level.

global gg
global dd

MutatedSCAgeFission = zeros(gg.numDiv,1);

% we start with all cells/mtDNA mutation free
MutatedAll = zeros(gg.numDiv,gg.initS);

% Find out the time at which the crypt fission event arose
CryptFisTime = size(s3);
CryptFisTime2 = CryptFisTime(1);

% Insert the data about the old crypt into the
MutatedAll(1:CryptFisTime(1),1:gg.initS) = s3;

% initiate time, time+1 means divTime has passed
time = CryptFisTime2;

%% Pre-determined random numbers for crypt simulation

% Mutation Rate random numbers
aaaa = DiscSampVec3((0:1),[gg.mutationRate1],(gg.mtDNA*gg.initS));

% Stem cell division type random numbers
bbbb = DiscSampVec2((1:3),[gg.Pa,gg.Ps,gg.Ps],gg.numDiv*gg.initS*2);
count2 = 1;

% Stem cell division type with advantage random numbers
cccc = DiscSampVec2((1:3),[gg.Pa-((gg.Ps*gg.adv)-gg.Ps),gg.Ps*gg.adv,gg.Ps],gg.numDiv*gg.initS*2);
count3 = 1;

% Segregation event random numbers
dddd = DiscSampVec2((1:2),[0.5,0.5],gg.numDiv*gg.initS*2);
count4 = 1;

% Stem cell replacement random numbers
eeee = rand(1,(gg.numDiv*gg.initS*2));
count5 = 1;

% Crypt Fission Crypt Numbering
ffff = rand(1,2*gg.numDiv);
count6 = 1;

%% Simulate only for certain time
while time < gg.numDiv
    
    if time == 1
        b = 0;
    else
        a = sum(MutatedAll(time,:));
        b = a > 0;
    end
    
    %% mutations occuring
    % random numbers generated for each mtDNA molecule
    % within all stem cells of the crypt to determine how many are
    % mutated
    
    if b == 0
        
        Mutated = [];
        
        for iii = 1:gg.initS
            
            Mutated(iii) = sum(aaaa(time,((iii*gg.mtDNA)-(gg.mtDNA-1)):...
                    ((iii*gg.mtDNA)-(gg.mtDNA-1)) + (gg.mtDNA-1)));
            
        end
        
        MutatedAll(time+1,:) = MutatedAll(time,:) + Mutated;
        
        time = time + 1;
        
    end
    
    if b > 0
        
        %% mutations occuring
        % random numbers generated for each mtDNA molecule
        % within all stem cells of the crypt to determine how many are
        % mutated
        
        Mutated = [];
        
        for iii = 1:gg.initS
            
           Mutated(iii) = sum(aaaa(time,((iii*gg.mtDNA)-(gg.mtDNA-1)):...
                    ((iii*gg.mtDNA)-(gg.mtDNA-1)) + (gg.mtDNA-1)));
                
        end
        
        MutatedAll(time+1,:) = MutatedAll(time,:) + Mutated;
        
        % if there are some mutations in our system, then we see how they
        % propagate
        RelevantMutations = find(MutatedAll(time,:)>0);
        
        % if there are mutations present
        if sum(MutatedAll(time,:))>0
            
            % for each cell with a mutation present
            for jj = 1 : numel(RelevantMutations)
                
                %%stem cell dividing (1 - asymmetric, 2 - symmetric 2 stem cells, 3 - symmetric 2 TA cells)
                
                if MutatedAll(time,RelevantMutations(jj)) > gg.mutThreshold*gg.mtDNA
                    divisionType = bbbb(count2);
                    count2 = count2 + 1;
                else
                    divisionType = cccc(count3);
                    count3 = count3 + 1;
                end
                
                % mutated mtDNA loss and gain before stem cell division
                % mutatedRep - how many new mutated mtDNAs you get in the stem
                % cell after doubling the number of mtDNA molecules.
                
                mutatedRep = DiscSampVec2...
                    ((0:gg.mtDNA),gg.RepProb...
                    (MutatedAll(time,RelevantMutations(jj)),:),1);
                
                % add new mtDNA mutation to old ones
                
                numMutated = mutatedRep + MutatedAll(time,RelevantMutations(jj));
                
                % division into two cells, each with n mtDNA
                % mutatedDiv - how many of the mutations will one cell get (the
                % other by proxy gets all the rest)
                
                mutatedDiv = DiscSampVec2...
                    ((0:gg.mtDNA),gg.DivProb...
                    (numMutated,:),1);
                
                % how many does the other cell have
                
                vectorDiv = [mutatedDiv, numMutated - mutatedDiv];
                
                % depending on the type of division, cells get kept or lost
                % asymmetric division occurs, one cell gets lost, one remains
                
                if divisionType == 1
                    remainingCell = dddd(count4);
                    count4 = count4 + 1;
                    remained = vectorDiv(remainingCell);
                    MutatedAll(time+1,RelevantMutations(jj)) = remained;
                    
                    % symmetric division into 2 stem cells, both are kept
                elseif divisionType == 2
                    remained1 = vectorDiv(1);
                    remained2 = vectorDiv(2);
                    MutatedAll(time+1,RelevantMutations(jj)) = remained1;
                    
                    % which one of the other cells will it replace?
                    a = 1:gg.initS;
                    possibleReplacements = a(a ~= RelevantMutations(jj));
                    b = ceil((gg.initS-1)*eeee(count5));
                    count5 = count5+1;
                    c = possibleReplacements(b);
                    MutatedAll(time+1,c) = remained2;
                    
                    % symmetric division into 2 TA cells, none are kept
                elseif divisionType == 3
                    
                    % which of the other ones gets doubled?
                    a = 1:gg.initS;
                    possibleReplacements = a(a ~= RelevantMutations(jj));
                    b = ceil((gg.initS-1)*eeee(count5));
                    count5 = count5+1;
                    c = possibleReplacements(b);
                    MutatedAll(time+1,RelevantMutations(jj)) = MutatedAll(time,c);
                end
            end
        end
        time = time + 1;
    end
end

%% Main Output

% How many stem cells at each age have a pathogenic mutation present.

for uu = 1 : gg.numDiv
    Mut = find(MutatedAll(uu,:) > (gg.mtDNA*gg.mutThreshold));
    MutNo = numel(Mut);
    MutatedSCAgeFission(uu,1) = MutNo;
end

%% Crypt Fission Events?

% The number of stem cells that are mutated in this crypt during its
% lifetime

SCMutatedNo = MutatedSCAgeFission(:,1);

% Primed scalar vector to record when and where a crypt fission event
% occurs

CryptFissionEvent = zeros(gg.numDiv,1);

% For each division that has occured, determine what the crypt
% fission probability is dependent on the number of stem cells that
% are mutated. Determine if fission does occur and record it in the
% CryptFissionEvent vector.

for hh = CryptFisTime2 : gg.numDiv
    
    SCMut = SCMutatedNo(hh,1);
    
    if SCMut == 0 && ffff(count6) < gg.cryptFissionProb
        CryptFissionEvent(hh,1) = 1;
    end
    
    count6 = count6 + 1;
    
    if SCMut > 0 &&  ffff(count6) < gg.cryptFissionProb*gg.cryptFissionFactor*SCMut
        CryptFissionEvent(hh,1) = 1;
    end
    
    count6 = count6 + 1;
    
end

if sum(CryptFissionEvent) > 0
    
    FissionAge = find(CryptFissionEvent == 1);
    
    for rr = 1 : numel(FissionAge)
        
        MutatedAllData = MutatedAll(1:FissionAge(rr),:);
        cryptFisSaveNo = ['dd.MutatedAllData' num2str(gg.cryptFisSave)];
        str = [cryptFisSaveNo, '= MutatedAllData;'];
        eval(str)
        gg.cryptFisSave = gg.cryptFisSave + 1;
    end
    
end

end


