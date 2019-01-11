function [MutatedSCAge, MutatedSCAgeCorr, MutatedSCAgeCorr2] = mtDNACrypt_ConstantV11FCN_COXAd()

% mtDNACrypt Function for Constant and Increasing Mutation Rate - FINAL

% The script is an amalgamation of the previous crypt model. It identifies
% that there are a certain number of mtDNA molecules residing within each
% stem cell of the crypt. With the evolution of stem cell divisions, the
% number of mutated mtDNA molecules evolves stochastically according to
% pre-determined probabilities. Also, with each additional mutated mtDNA
% molecule, the model determines which kind of mutation has developed
% according to probability data previously aquired. Therefore, this model
% is a more accurate representation of the processes that take place within
% the crypt and at the tissue level.

% Tracks multiple mutations on single mtDNA species for clonal expansion
% comparison of multiple mutations within individual cells with biological
% data

% Bring in the global variable 'gg' that has already been set up and make a
% global 'dd' variable
global gg
global dd

% Set up the results matrices
MutatedSCAge = zeros(gg.numDiv,gg.numRuns);
MutatedSCAgeCorr = zeros(gg.numDiv,gg.numRuns);
MutatedSCAgeCorr2 = zeros(gg.numDiv,gg.numRuns);

% Set up a progress tracking bar
h = waitbar(0,'Simulating model w/o crypt fission, please wait...');

% Set up the multiple mutations record matrices for stem cells at age 70
% years of simulated time
SingleMutRecord = zeros(gg.numRuns,gg.initS);
MultipleMutRecord = zeros(gg.numRuns,gg.initS);

% Probability for each age (numDiv) getting a mutation
mutProbAge = zeros(2,gg.numDiv);

for pp = 1 : gg.numRuns
    
    % set up the multiple mutations species ID result structure
    mtDNAmutations = zeros(gg.numDiv, gg.initS*gg.mtDNA);
    
    % we start with all cells/mtDNA mutation free
    MutatedAll = zeros(gg.numDiv,gg.initS);
    
    % Set up the first value of the original mutation
    origMut = 1;
    
    % Set up the mtDNA species records
    speciesIDRecord = [];
    speciesIDMultRecord = [];
    
    % initiate time, time+1 means divTime has passed
    time = 1;
    
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
    
    % Species ID checking
    ffff = ceil(rand(1,(gg.mtDNA*gg.numDiv))*gg.mtDNA);
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
            
            % For each number of new mutations, determine if it is
            % replacing any of the current mutated mtDNA molecules. If it
            % is replacing any, "Mutated" is decreased by the same amount
            % for that stem cell.
            
            % Proceed if there are mutations present
            
            MutatedAll(time+1,:) = MutatedAll(time,:) + Mutated;
            
            % This is the point at which the first mutation will emerge
            % First mutation needs to be inserted and recorded
            % This is just mutation insertion only where they appear
            
            if max(Mutated) > 0
                
                for iii = 1 : gg.initS
                    
                    % Records how many mtDNA acquire second mutation per
                    % stem cell
                    
                    Multiple = 0;
                    
                    if Mutated(iii) > 0
                        
                        % Isolate the current stem cells mutational species
                        
                        tempA = mtDNAmutations(time, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA)));
                        
                        % Find all the WT mtDNA molecules
                        
                        tempC = find(tempA == 0);
                        
                        % For each new mutation, determine whether it is
                        % affecting a WT mtDNA or an already mutated mtDNA
                        % molecule
                        
                        ttt = 1;
                        mutPos = zeros(1,Mutated(iii));
                        
                        while ttt <= Mutated(iii)
                            
                            tempZ = ffff(count6);
                            
                            if isempty(find(mutPos == tempZ))
                                
                                mutPos(ttt) = tempZ;
                                count6 = count6 + 1;
                                ttt = ttt + 1;
                                
                            else
                                
                                count6 = count6 + 1;
                                
                            end
                            
                        end
                        
                        % Which values of mutPos are not present in tempC
                        
                        for ttt = 1 : numel(mutPos)
                            occuMut = find(tempC == mutPos(ttt));
                            if isempty(occuMut)
                                Multiple = Multiple + 1;
                            end
                        end
                        
                        if Multiple > 0
                            
                            % Find all current mutations
                            
                            currMut = find(tempA > 0);
                            
                            % Produce a random permutation of the indexed
                            % mutated mtDNA molecules
                            
                            currMutRandom = currMut(randperm(numel(currMut)));
                            
                            % The overwritten molecules species ID's will be
                            
                            overSpeciesID = tempA(currMutRandom(1:Multiple));
                            
                            % Determine the new species IDs for the mutated mtDNA molecules
                            
                            tempE = origMut : (origMut + (Multiple-1));
                            
                            % Take away the multiple mutations from the species
                            % ID generator vector
                            
                            tempEMult = tempE(1 : Multiple); % Contains the multiple species ID
                            
                            tempEWT = tempE((Multiple+1) : end); % Contains the normal species ID
                            
                            % Insert the new species ID into the current list
                            % of species IDs in the WT molecule positions
                            
                            tempA(tempC(1:numel(tempEWT))) = tempEWT;
                            
                            % Replace the selected mtDNA species to be
                            % overwritten for multiple mutations on same
                            % mtDNA species.
                            
                            for ttt =1 : Multiple
                                a = find(tempA == overSpeciesID(ttt));
                                tempA(a(end)) = tempEMult(ttt);
                            end
                            
                            % Insert the modified species ID vector into the
                            % master matrix
                            
                            mtDNAmutations(time+1, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA))) = tempA;
                            
                            % This needs to be reflected in the MutatedAll
                            % array as well
                            
                            MutatedAll(time+1,iii) = MutatedAll(time+1,iii) - Multiple;
                            
                            % increase the species ID tracker by one
                            
                            origMut = origMut + numel(tempE);
                            
                            % Recording the multiple mutation information
                            % For each mutation, need to check whether it
                            % has been muliplied before.
                            
                            for ttt = 1 : Multiple
                                
                                % determine whether the speciesID that is about
                                % to be overwritten has any multiple mutations
                                % already
                                
                                a = find(speciesIDRecord == overSpeciesID(ttt));
                                
                                if isempty(a)
                                    speciesIDRecord(end+1) = tempEMult(ttt);
                                    speciesIDMultRecord(end+1) = 2;
                                else
                                    speciesIDRecord(end+1) = tempEMult(ttt);
                                    speciesIDMultRecord(end+1) = speciesIDMultRecord(a) + 1;
                                end
                            end
                            
                        else
                            
                            % Determine the new species IDs for the mutated mtDNA molecules
                            
                            tempE = origMut : (origMut + Mutated(iii)-1);
                            
                            % Insert the new species ID into the current list
                            % of species IDs in the WT molecule positions
                            
                            tempA(tempC(1:numel(tempE))) = tempE;
                            
                            % Insert the modified species ID vector into the
                            % master matrix
                            
                            mtDNAmutations(time+1, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA))) = tempA;
                            
                            % increase the species ID tracker by one
                            
                            origMut = origMut + numel(tempE);
                            
                        end
                        
                    end
                    
                end
                
            end
            
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
            
            MutatedAll(time,:) = MutatedAll(time,:) + Mutated;
            
            % This is the point at which additional mutations will arise
            % and where multiple mutations will be tracked and recorded
            
            if max(Mutated) > 0
                
                for iii = 1 : gg.initS
                    
                    % Records how many mtDNA acquire second mutation per
                    % stem cell
                    
                    Multiple = 0;
                    
                    if Mutated(iii) > 0
                        
                        % Isolate the current stem cells mutational species
                        
                        tempA = mtDNAmutations(time, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA)));
                        
                        % Find all the WT mtDNA molecules
                        
                        tempC = find(tempA == 0);
                        
                        % For each new mutation, determine whether it is
                        % affecting a WT mtDNA or an already mutated mtDNA
                        % molecule
                        
                        ttt = 1;
                        mutPos = zeros(1,Mutated(iii));
                        
                        while ttt <= Mutated(iii)
                            
                            tempZ = ffff(count6);
                            
                            if isempty(find(mutPos == tempZ)) % For same number sequence in ffff check
                                
                                mutPos(ttt) = tempZ;
                                count6 = count6 + 1;
                                ttt = ttt + 1;
                                
                            else
                                
                                count6 = count6 + 1;
                                
                            end
                            
                        end
                        
                        % Which values of mutPos are not present in tempC
                        
                        for ttt = 1 : numel(mutPos)
                            occuMut = find(tempC == mutPos(ttt));
                            if isempty(occuMut)
                                Multiple = Multiple + 1;
                            end
                        end
                        
                        if Multiple > 0
                            
                            % Find all current mutations
                            
                            currMut = find(tempA > 0);
                            
                            % Produce a random permutation of the indexed
                            % mutated mtDNA molecules
                            
                            currMutRandom = currMut(randperm(numel(currMut)));
                            
                            % The overwritten molecules species ID's will be
                            
                            overSpeciesID = tempA(currMutRandom(1:Multiple));
                            
                            % Determine the new species IDs for the mutated mtDNA molecules
                            
                            tempE = origMut : (origMut + (Multiple-1));
                            
                            % Take away the multiple mutations from the species
                            % ID generator vector
                            
                            tempEMult = tempE(1 : Multiple); % Contains the multiple species ID
                            
                            tempEWT = tempE((Multiple+1) : end); % Contains the normal species ID
                            
                            % Insert the new species ID into the current list
                            % of species IDs in the WT molecule positions
                            
                            tempA(tempC(1:numel(tempEWT))) = tempEWT;
                            
                            % Replace the selected mtDNA species to be
                            % overwritten for multiple mutations on same
                            % mtDNA species.
                            
                            for ttt = 1 : Multiple
                                a = find(tempA == overSpeciesID(ttt));
                                tempA(a(end)) = tempEMult(ttt);
                            end
                            
                            % insert the modified species ID vector into the
                            % master matrix
                            
                            mtDNAmutations(time, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA))) = tempA;
                            
                            % This needs to be reflected in the MutatedAll
                            % Array as well
                            
                            MutatedAll(time,iii) = MutatedAll(time,iii) - Multiple;
                            
                            % Increase the species ID tracker by one
                            
                            origMut = origMut + numel(tempE);
                            
                            % Recording the multiple mutation information
                            % For each mutation, need to check whether it
                            % has been muliplied before.
                            
                            for ttt = 1 : Multiple
                                % determine whether the speciesID that is about
                                % to be overwritten has any multiple mutations
                                % already
                                
                                a = find(speciesIDRecord == overSpeciesID(ttt));
                                
                                if isempty(a)
                                    speciesIDRecord(end+1) = tempEMult(ttt);
                                    speciesIDMultRecord(end+1) = 2;
                                else
                                    speciesIDRecord(end+1) = tempEMult(ttt);
                                    speciesIDMultRecord(end+1) = speciesIDMultRecord(a) + 1;
                                end
                            end
                            
                        else
                            
                            % Determine the new species IDs for the mutated mtDNA molecules
                            
                            tempE = origMut : (origMut + Mutated(iii)-1);
                            
                            % Insert the new species ID into the current list
                            % of species IDs in the WT molecule positions
                            
                            tempA(tempC(1:numel(tempE))) = tempE;
                            
                            % insert the modified species ID vector into the
                            % master matrix
                            
                            mtDNAmutations(time, (((iii*gg.mtDNA)-(gg.mtDNA-1)):(iii*gg.mtDNA))) = tempA;
                            
                            % increase the species ID tracker by one
                            
                            origMut = origMut + numel(tempE);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            % if there are some mutations in our system, then we see how they
            % propagate
            
            RelevantMutations = find(MutatedAll(time,:)>0);
            
            % if there are mutations present
            
            if sum(MutatedAll(time,:))>0
                
                % for each cell with a mutation present
                
                for jj = 1 : numel(RelevantMutations)
                    
                    % stem cell dividing (1 - asymmetric, 2 - symmetric 2 stem cells, 3 - symmetric 2 TA cells)
                    
                    if MutatedAll(time,RelevantMutations(jj)) >= gg.mutThreshold*gg.mtDNA
                        divisionType = bbbb(count2);
                        count2 = count2 + 1;
                    else
                        divisionType = cccc(count3);
                        count3 = count3 + 1;
                    end
                    
                    % mutated mtDNA loss and gain before stem cell division
                    % mutratedRep - how many new mutated mtDNAs you get in the stem
                    % cell after doubling the number of mtDNA molecules.
                    
                    mutatedRep = DiscSampVec2...
                        ((0:gg.mtDNA),gg.RepProb...
                        (MutatedAll(time,RelevantMutations(jj)),:),1);
                    
                    % add new mtDNA mutation to old ones
                    
                    numMutated = mutatedRep + MutatedAll(time,RelevantMutations(jj));
                    
                    % At this point the multiple mutations in
                    % mtDNAmutations need to be increased to
                    % the numbers that are in numMutated
                    
                    tempA = numMutated;
                    tempB = mtDNAmutations(time,...
                        (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                        (RelevantMutations(jj)*gg.mtDNA)));
                    tempC = tempB(tempB>0);
                    tempD = tempC(randi(numel(tempC),1,tempA));
                    
                    % Store this matrix to a seperate variable
                    
                    numMutatedMutations = tempD;
                    
                    % division into two cells, each with n mtDNA
                    % mutatedDiv - how many of the mutations will one cell
                    % get (the other by proxy gets all the rest)
                    
                    % Altered for DivProb with advantage...
                    
                    if divisionType == 1
                        
                        mutatedDiv = DiscSampVec2...
                            ((0:gg.mtDNA),gg.DivProb103...
                            (numMutated,:),1);
                        
                    else
                        
                        mutatedDiv = DiscSampVec2...
                            ((0:gg.mtDNA),gg.DivProb...
                            (numMutated,:),1);
                        
                    end
                    
                    % Depeding on the number of mtDNA molecules go into one
                    % cell, the other gets the other lot this is based on
                    % numMutatedMutations in the master cell before
                    % segregation
                    
                    tempA = mutatedDiv;
                    tempB = numMutatedMutations(randperm(numel(numMutatedMutations)));
                    
                    Cell1 = tempB(1:tempA);
                    Cell2 = tempB(tempA+1 : end);
                    
                    if isempty(Cell1)
                        Cell1 = 0;
                    end
                    
                    if isempty(Cell2)
                        Cell2 = 0;
                    end
                    
                    % how many does the other cell have
                    
                    vectorDiv = [mutatedDiv, numMutated - mutatedDiv];
                    
                    % depending on the type of division, cells get kept or lost
                    % asymmetric division occurs, one cell gets lost, one remains
                    
                    if divisionType == 1
                        
                        remainingCell = 2; % Advantage forces the mutatedDiv result to be the stem cell
                        count4 = count4 + 1;
                        remained = vectorDiv(remainingCell);
                        MutatedAll(time+1,RelevantMutations(jj)) = remained;
                        
                        if remainingCell == 1
                            
                            % insert the new cell multiple mutation data
                            % depending on which cell is chosen for an
                            % aysmmetric division fate outcome
                            
                            tempA = zeros(1,gg.mtDNA);
                            tempA(1:numel(Cell1)) = Cell1;
                            mtDNAmutations(time+1,...
                                (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                            
                        else
                            
                            tempA = zeros(1,gg.mtDNA);
                            tempA(1:numel(Cell2)) = Cell2;
                            mtDNAmutations(time+1,...
                                (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                            
                        end
                        
                        % symmetric division into 2 stem cells, both are kept
                        
                    elseif divisionType == 2
                        
                        remained1 = vectorDiv(1);
                        remained2 = vectorDiv(2);
                        MutatedAll(time+1,RelevantMutations(jj)) = remained1;
                        
                        % insert the new cell multiple mutation data for
                        % the stem cell that stays for the symmetric fate
                        % outcome 1
                        
                        tempA = zeros(1,gg.mtDNA);
                        tempA(1:numel(Cell1)) = Cell1;
                        mtDNAmutations(time+1,...
                            (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                            (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                        
                        % which one of the other cells will it replace?
                        
                        a = 1:gg.initS;
                        possibleReplacements = a(a ~= RelevantMutations(jj));
                        b = ceil((gg.initS-1)*eeee(count5));
                        count5 = count5+1;
                        c = possibleReplacements(b);
                        MutatedAll(time+1,c) = remained2;
                        
                        % insert the new cell multiple mutation data for
                        % the stem cell that stays for the symmetric fate
                        % outcome 2
                        
                        tempA = zeros(1,gg.mtDNA);
                        tempA(1:numel(Cell2)) = Cell2;
                        mtDNAmutations(time+1,...
                            (((c*gg.mtDNA)-(gg.mtDNA-1)):...
                            (c*gg.mtDNA))) = tempA;
                        
                        % symmetric division into 2 TA cells, none are kept
                        
                    elseif divisionType == 3
                        
                        % which of the other ones gets doubled?
                        
                        a = 1:gg.initS;
                        possibleReplacements = a(a ~= RelevantMutations(jj));
                        b = ceil((gg.initS-1)*eeee(count5));
                        count5 = count5+1;
                        c = possibleReplacements(b);
                        MutatedAll(time+1,RelevantMutations(jj)) = MutatedAll(time,c);
                        
                        % insert the new cell multiple mutation data for
                        % the stem cell that stays for the symmetric fate
                        % outcome 2
                        
                        mtDNAmutations(time+1,(((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                            (RelevantMutations(jj)*gg.mtDNA))) = ...
                            mtDNAmutations(time,(((c*gg.mtDNA)-(gg.mtDNA-1)):...
                            (c*gg.mtDNA)));
                        
                    end
                    
                end
                
            end
            
            %%%%%%%%%%%% COX DEF SC RATE INCREASE %%%%%%%%%%%%%%%%%%%%
            
            if strcmp(gg.COXSCTimePoint, 'Yes') == 1
                
                % Run just the DIVISION CODE again for COX neg stem cells at
                % specific time points
                
                % Run code every n timepoints
                
                if mod(time,gg.COXSCTimePointInterval) == 0
                    
                    COXDefCycle = 0;
                    
                    while COXDefCycle < gg.COXDefCycleRepeats
                        
                        blueSCPres = find(MutatedAll(time+1,:)>=(gg.mtDNA*gg.mutThreshold));
                        
                        if ~isempty(blueSCPres)
                            
                            RelevantMutations = blueSCPres;
                            
                            % if there are mutations present
                            if sum(MutatedAll(time+1,:))>0
                                
                                % for each cell with a mutation present
                                for jj = 1 : numel(RelevantMutations)
                                    
                                    % stem cell dividing (1 - asymmetric, 2 - symmetric 2 stem cells, 3 - symmetric 2 TA cells)
                                    
                                    if MutatedAll(time+1,RelevantMutations(jj)) >= gg.mutThreshold*gg.mtDNA
                                        divisionType = bbbb(count2);
                                        count2 = count2 + 1;
                                    else
                                        divisionType = cccc(count3);
                                        count3 = count3 + 1;
                                    end
                                    
                                    % mutated mtDNA loss and gain before stem cell division
                                    % mutratedRep - how many new mutated mtDNAs you get in the stem
                                    % cell after doubling the number of mtDNA molecules.
                                    
                                    mutatedRep = DiscSampVec2...
                                        ((0:gg.mtDNA),gg.RepProb...
                                        (MutatedAll(time+1,RelevantMutations(jj)),:),1);
                                    
                                    % add new mtDNA mutation to old ones
                                    
                                    numMutated = mutatedRep + MutatedAll(time+1,RelevantMutations(jj));
                                    
                                    % At this point the multiple mutations in
                                    % mtDNAmutations need to be increased to
                                    % the numbers that are in numMutated
                                    
                                    tempA = numMutated;
                                    tempB = mtDNAmutations(time+1,...
                                        (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                        (RelevantMutations(jj)*gg.mtDNA)));
                                    tempC = tempB(tempB>0);
                                    tempD = tempC(randi(numel(tempC),1,tempA));
                                    
                                    % Store this matrix to a seperate variable
                                    
                                    numMutatedMutations = tempD;
                                    
                                    % division into two cells, each with n mtDNA
                                    % mutatedDiv - how many of the mutations will one cell get (the
                                    % other by proxy gets all the rest)
                                    
                                    % Altered for DivProb with advantage...
                                    
                                    if divisionType == 1
                                        
                                        mutatedDiv = DiscSampVec2...
                                            ((0:gg.mtDNA),gg.DivProb103...
                                            (numMutated,:),1);
                                        
                                    else
                                        
                                        mutatedDiv = DiscSampVec2...
                                            ((0:gg.mtDNA),gg.DivProb...
                                            (numMutated,:),1);
                                        
                                    end
                                    
                                    % Depeding on the number of mtDNA molecules go into one
                                    % cell, the other gets the other lot this is based on
                                    % numMutatedMutations in the master cell before
                                    % segregation
                                    
                                    tempA = mutatedDiv;
                                    tempB = numMutatedMutations(randperm(numel(numMutatedMutations)));
                                    
                                    Cell1 = tempB(1:tempA);
                                    Cell2 = tempB(tempA+1 : end);
                                    
                                    if isempty(Cell1)
                                        Cell1 = 0;
                                    end
                                    
                                    if isempty(Cell2)
                                        Cell2 = 0;
                                    end
                                    
                                    % how many does the other cell have
                                    
                                    vectorDiv = [mutatedDiv, numMutated - mutatedDiv];
                                    
                                    % depending on the type of division, cells get kept or lost
                                    % asymmetric division occurs, one cell gets lost, one remains
                                    
                                    if divisionType == 1
                                        
                                        remainingCell =  2; % Advantage forces the mutatedDiv result to be the stem cell
                                        count4 = count4 + 1;
                                        remained = vectorDiv(remainingCell);
                                        MutatedAll(time+1,RelevantMutations(jj)) = remained;
                                        
                                        if remainingCell == 1
                                            
                                            % insert the new cell multiple mutation data
                                            % depending on which cell is chosen for an
                                            % aysmmetric division fate outcome
                                            
                                            tempA = zeros(1,gg.mtDNA);
                                            tempA(1:numel(Cell1)) = Cell1;
                                            mtDNAmutations(time+1,...
                                                (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                                (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                                            
                                        else
                                            
                                            tempA = zeros(1,gg.mtDNA);
                                            tempA(1:numel(Cell2)) = Cell2;
                                            mtDNAmutations(time+1,...
                                                (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                                (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                                            
                                        end
                                        
                                        % symmetric division into 2 stem cells, both are kept
                                        
                                    elseif divisionType == 2
                                        
                                        remained1 = vectorDiv(1);
                                        remained2 = vectorDiv(2);
                                        MutatedAll(time+1,RelevantMutations(jj)) = remained1;
                                        
                                        % insert the new cell multiple mutation data for
                                        % the stem cell that stays for the symmetric fate
                                        % outcome 1
                                        
                                        
                                        tempA = zeros(1,gg.mtDNA);
                                        tempA(1:numel(Cell1)) = Cell1;
                                        mtDNAmutations(time+1,...
                                            (((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                            (RelevantMutations(jj)*gg.mtDNA))) = tempA;
                                        
                                        % which one of the other cells will it replace?
                                        
                                        a = 1:gg.initS;
                                        possibleReplacements = a(a ~= RelevantMutations(jj));
                                        b = ceil((gg.initS-1)*eeee(count5));
                                        count5 = count5+1;
                                        c = possibleReplacements(b);
                                        MutatedAll(time+1,c) = remained2;
                                        
                                        % insert the new cell multiple mutation data for
                                        % the stem cell that stays for the symmetric fate
                                        % outcome 2
                                        
                                        tempA = zeros(1,gg.mtDNA);
                                        tempA(1:numel(Cell2)) = Cell2;
                                        mtDNAmutations(time+1,...
                                            (((c*gg.mtDNA)-(gg.mtDNA-1)):...
                                            (c*gg.mtDNA))) = tempA;
                                        
                                        % symmetric division into 2 TA cells, none are kept
                                        
                                    elseif divisionType == 3
                                        
                                        % which of the other ones gets doubled?
                                        
                                        a = 1:gg.initS;
                                        possibleReplacements = a(a ~= RelevantMutations(jj));
                                        b = ceil((gg.initS-1)*eeee(count5));
                                        count5 = count5+1;
                                        c = possibleReplacements(b);
                                        MutatedAll(time+1,RelevantMutations(jj)) = MutatedAll(time+1,c);
                                        
                                        % insert the new cell multiple mutation data for
                                        % the stem cell that stays for the symmetric fate
                                        % outcome 2
                                        
                                        mtDNAmutations(time+1,(((RelevantMutations(jj)*gg.mtDNA)-(gg.mtDNA-1)):...
                                            (RelevantMutations(jj)*gg.mtDNA))) = ...
                                            mtDNAmutations(time+1,(((c*gg.mtDNA)-(gg.mtDNA-1)):...
                                            (c*gg.mtDNA)));
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                        end
                        
                        COXDefCycle = COXDefCycle + 1;
                        
                    end
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            time = time + 1;
            
        end
        
    end
    
    %% METRICS mtDNA clonal expansion
    
    % Identifying successul and failed mtDNA clonal expansion events.
    
    % Insert an extra column into MutatedAll so if a mutation appears at the
    % first time point then the difference is captured
    
    MutBuffer = zeros(1,gg.initS);
    
    MutatedAll2 = [MutBuffer; MutatedAll];
    
    % Find the difference between mutation events and clonal expansion
    
    Difference = zeros(gg.numDiv+1,gg.initS);
    
    for ii = 1 : gg.initS % InitS
        for tt = 1 : gg.numDiv-1
            Difference(tt+1,ii) = MutatedAll2(tt+1,ii) - MutatedAll2(tt,ii);
        end
    end
    
    % Find where the 1 mtDNA values are for all stem cells in the niche
    
    for ll = 1 : gg.initS
        
        PrimaryName = ['PrimaryMut' num2str(ll)];
        
        PrimaryMut = find(MutatedAll2(:,ll) > 0)';
        
        str = [PrimaryName, '=PrimaryMut;'];
        
        eval(str)
        
        a = eval(['PrimaryMut' num2str(ll)]);
        
        for uu = 1 : numel(a);
            
            if Difference(a(uu),ll) == MutatedAll2(a(uu),ll);
                a(uu) = a(uu);
            else
                a(uu) = 0;
            end
            
        end
        
        a(a == 0) = [];
        
        str = [PrimaryName, '=a;'];
        eval(str)
        
    end
    
    % Find where the end of the clonal expansion is if it does have an end
    
    for ll = 1 : gg.initS
        
        PrimaryName = ['PrimaryEndMut' num2str(ll)];
        
        PrimaryEndMut = find(MutatedAll2(:,ll) == 0)';
        
        str = [PrimaryName, '=PrimaryEndMut;'];
        
        eval(str)
        
        a = eval(['PrimaryEndMut' num2str(ll)]);
        
        for uu = 2 : numel(a);
            
            if Difference(a(uu),ll) == MutatedAll2(a(uu)-1,ll)*-1 && Difference(a(uu),ll) ~= 0
                a(uu) = a(uu);
            else
                a(uu) = 0;
            end
            
        end
        
        a(a == 0) = [];
        a(a == 1) = [];
        
        str = [PrimaryName, '=a;'];
        eval(str)
        
    end
    
    % Find the failed clonal expansions and the times
    
    for ll = 1:gg.initS
        
        a = eval(['PrimaryMut' num2str(ll)]);
        b = eval(['PrimaryEndMut' num2str(ll)]);
        
        start = numel(a);
        finish = numel(b);
        
        if start == finish
            for tt = 1 : start
                vv = b(tt) - a(tt);
                if vv >= 0
                    gg.FailedCE(1,end+1) = vv;
                    gg.FailedCE(2,end) = a(tt);
                    gg.FailedCE(3,end) = b(tt);
                end
            end
        end
        
        if start > finish
            
            for jj = 1 : finish
                gg.FailedCE(1,end+1) = b(jj) - a(jj);
                gg.FailedCE(2,end) = a(jj);
                gg.FailedCE(3,end) = b(jj);
            end
            
            for jj = start
                c = find(MutatedAll2(a(jj):gg.numDiv,ll) == gg.mtDNA);
                d = min(c) + a(jj) - 1;
                
                if isempty(c) % For those that are still transient
                else
                    vv = d - a(jj);
                    if vv >= 0
                        gg.SuccessCE(1,end+1) = vv + 1;
                        gg.SuccessCE(2,end) = a(jj);
                        gg.SuccessCE(3,end) = d;
                    end
                end
                
            end
        end
    end
    
    %% METRICS SC Niche Succession
    
    stemCellAll = zeros(gg.numDiv+1,1);
    
    for ii = 1 : gg.numDiv
        stemCellAll(ii+1,1) = numel(find(MutatedAll(ii,:) >= gg.mutThreshold*gg.mtDNA));
    end
    
    % Find the difference between mutation SC and niche succession
    
    SCDifference = zeros(gg.numDiv+1,1);
    
    for tt = 1 : gg.numDiv - 1
        SCDifference(tt+1,1) = stemCellAll(tt+1,1) - stemCellAll(tt,1);
    end
    
    % Make both the Difference and the MutatedAll the same size
    
    SCPrimaryMut = find(stemCellAll > 0)';
    
    for uu = 1 : numel(SCPrimaryMut);
        
        if SCDifference(SCPrimaryMut(uu),1) == stemCellAll(SCPrimaryMut(uu),1);
            SCPrimaryMut(uu) = SCPrimaryMut(uu);
        else
            SCPrimaryMut(uu) = 0;
        end
        
    end
    
    SCPrimaryMut(SCPrimaryMut == 0) = [];
    
    % Find where the end of the clonal expansion is if it does have an end
    
    SCPrimaryEndMut = find(stemCellAll == 0)';
    
    for uu = 2 : numel(SCPrimaryEndMut);
        
        if SCDifference(SCPrimaryEndMut(uu),1) == stemCellAll(SCPrimaryEndMut(uu)-1,1)*-1 && SCDifference(SCPrimaryEndMut(uu),1) ~= 0
            SCPrimaryEndMut(uu) = SCPrimaryEndMut(uu);
        else
            SCPrimaryEndMut(uu) = 0;
        end
        
    end
    
    SCPrimaryEndMut(SCPrimaryEndMut == 0) = [];
    SCPrimaryEndMut(SCPrimaryEndMut == 1) = [];
    
    
    % Find the failed clonal expansions and the times
    
    start = numel(SCPrimaryMut);
    finish = numel(SCPrimaryEndMut);
    
    if start == finish
        for tt = 1 : start
            gg.NicheFailedSC(1,end+1) = SCPrimaryEndMut(tt) - SCPrimaryMut(tt);
            gg.NicheFailedSC(2,end) = SCPrimaryMut(tt);
            gg.NicheFailedSC(3,end) = SCPrimaryEndMut(tt);
        end
    end
    
    if start > finish
        for jj = 1 : finish
            gg.NicheFailedSC(1,end+1) = SCPrimaryEndMut(jj) - SCPrimaryMut(jj);
            gg.NicheFailedSC(2,end) = SCPrimaryMut(jj);
            gg.NicheFailedSC(3,end) = SCPrimaryEndMut(jj);
        end
        for jj = start
            c = find(stemCellAll(SCPrimaryMut(jj):gg.numDiv,1) == gg.initS);
            d = min(c) + SCPrimaryMut(jj) - 1;
            
            if isempty(c)
            else
                vv = d - SCPrimaryMut(jj);
                if vv >= 0
                    gg.NicheSuccessSC(1,end+1) = d - SCPrimaryMut(jj) + 1;
                    gg.NicheSuccessSC(2,end) = SCPrimaryMut(jj) + 1;
                    gg.NicheSuccessSC(3,end) = d;
                end
            end
        end
    end
    
    %% Correction factor for mtDNAmutations and mutatedAll
    
    % This part of the code affects both mtDNAmutations and mutatedAll in
    % order to affect MutatedSCAge to determine how much COX deficiency
    % will be present after the correction factor has been implemented.
    % This will run alongside the current code so there is a measure of the
    % affect the correction factor has
    
    % Determine the max number of mutations present
    
    maxMut = max(max(mtDNAmutations));
    
    % Generate each speciesID
    
    maxSpecies = 1 : maxMut;
    
    % For every species ID thats present in specesIDRecord, delete from
    % maxSpecies
    
    for vv = 1 : numel(speciesIDRecord)
        maxSpecies(maxSpecies == speciesIDRecord(vv)) = [];
    end
    
    % Determine which numbers need to be excluded from the list present,
    % need to use a random number generator
    
    maxRand = rand(1,numel(maxSpecies));
    
    corrPos = maxRand <= gg.COXCorrectionFactor;
    
    exclSpecies1 = maxSpecies(corrPos);
    
    % Now for the species that have multiple mutations present. Each
    % mutation has to be assessed individually
    
    % Generate the number of random numbers required for each mutation
    
    multRand = rand(1,sum(speciesIDMultRecord));
    
    % Go through each species with multiple mutations and see if any
    % dont contain any COX deficiency mutation
    
    exclSpecies2 = [];
    
    for vv = 1 : numel(speciesIDRecord)
        
        if min(multRand(1:speciesIDMultRecord(vv))) <= (gg.COXCorrectionFactor)
            
            exclSpecies2(end+1) = speciesIDRecord(vv);
            
        end
        
        multRand(1:speciesIDMultRecord(vv)) = [];
        
    end
    
    % Combine both exclSpecies and exclSpecies2 which contain the species
    % IDs that are to be excluded from mtDNAmutations.
    
    exclSpecies = [exclSpecies1 exclSpecies2];
    
    % Delete the numbers that are present in corrPos from mtDNAmutations
    
    mtDNAmutationsCorr = mtDNAmutations;
    
    for vv = 1 : numel(exclSpecies)
        mtDNAmutationsCorr(mtDNAmutationsCorr == exclSpecies(vv)) = 0;
    end
    
    %% Correction factor for mtDNAmutations and mutatedAll - Adjusted
    
    % This part of the code affects both mtDNAmutations and mutatedAll in
    % order to affect MutatedSCAge to determine how much COX deficiency
    % will be present after the correction factor has been implemented.
    % This will run alongside the current code so there is a measure of the
    % affect the correction factor has.
    
    % Set up the vector that is going to record the mtDNAspecies that are
    % homoplasmic within the cell.
    
    homoplas_mtDNASpecies = [];
    
    % For each age and for each stem cell determine where the homoplasmic
    % mutations are.
    
    for vv = 1 : gg.numDiv
        
        for bb = 1 : gg.initS
            
            % Determine the vector to be assessed (stem cell at timepoint)
            
            vectorCorr = mtDNAmutationsCorr(vv, (((bb*gg.mtDNA)-(gg.mtDNA-1)):(bb*gg.mtDNA)));
            
            % What are the unique values present within this
            
            vectorCorr_unique = unique(vectorCorr);
            
            vectorCorr_unique(vectorCorr_unique == 0) = [];
            
            % For each unique speciesID, what is the %
            
            for jj = 1 : numel(vectorCorr_unique)
                
                vectorCorr_number = numel(find(vectorCorr == vectorCorr_unique(jj)));
                vectorCorr_percentage = vectorCorr_number / gg.mtDNA * 100;
                
                if vectorCorr_percentage == 100
                    
                    % Set what happens when there is a homoplasmic mtDNA
                    % species present -- It gets recorded into a new vector
                    
                    homoplas_mtDNASpecies(end+1) = vectorCorr_unique(jj);
                    
                end
                
            end
            
        end
        
    end
    
    % Need to get rid of repeated values in order
    
    homoplas_mtDNASpecies = unique(homoplas_mtDNASpecies);
    
    %% Need to remove the species IDs that dont satisfy the inclusion criteria
    
    homoplas_mtDNASpecies_post = homoplas_mtDNASpecies(rand(1,...
        numel(homoplas_mtDNASpecies)) <= gg.COXCorrectionFactor2);
    
    % Delete the numbers that are present in homoplas_mtDNASpecies_post from mtDNAmutations
    
    mtDNAmutationsCorr2 = mtDNAmutationsCorr;
    
    for vv = 1 : numel(homoplas_mtDNASpecies_post)
        mtDNAmutationsCorr2(mtDNAmutationsCorr2 == homoplas_mtDNASpecies_post(vv)) = 0;
    end
    
    
    %% Correction Factor Integration
    
    % Now that the correction factor has been implemented, we need to
    % determine, for each age, for each stem cell, the new number of mtDNA
    
    % mtDNAmutationsCorr summed up in MutatedAllCorr
    
    MutatedAllCorr = zeros(gg.numDiv,gg.initS);
    
    for vv = 1 : gg.numDiv
        
        for uu = 1 : gg.initS
            
            section = mtDNAmutationsCorr(vv,((gg.mtDNA*uu) - (gg.mtDNA-1)) : (gg.mtDNA*uu));
            speciesPres = find(section > 0);
            numSpeciesPresent = numel(speciesPres);
            MutatedAllCorr(vv,uu) = numSpeciesPresent;
            
        end
        
    end
    
    % Main Output for correctionFactorResult
    
    for uu = 1 : gg.numDiv
        Mut = find(MutatedAllCorr(uu,:) >= (gg.mtDNA*gg.mutThreshold));
        MutNo = numel(Mut);
        MutatedSCAgeCorr(uu,pp) = MutNo;
    end
    
    %% Correction Factor 2 Integration
    
    % Now that the correction factor has been implemented, we need to
    % determine, for each age, for each stem cell, the new number of mtDNA
    
    % mtDNAmutationsCorr2 summed up in MutatedAllCorr2
    
    MutatedAllCorr2 = zeros(gg.numDiv,gg.initS);
    
    for vv = 1 : gg.numDiv
        
        for uu = 1 : gg.initS
            
            section = mtDNAmutationsCorr2(vv,((gg.mtDNA*uu) - (gg.mtDNA-1)) : (gg.mtDNA*uu));
            speciesPres = find(section > 0);
            numSpeciesPresent = numel(speciesPres);
            MutatedAllCorr2(vv,uu) = numSpeciesPresent;
            
        end
        
    end
    
    % Main Output for correctionFactorResult
    
    for uu = 1 : gg.numDiv
        Mut = find(MutatedAllCorr2(uu,:) >= (gg.mtDNA*gg.mutThreshold));
        MutNo = numel(Mut);
        MutatedSCAgeCorr2(uu,pp) = MutNo;
    end
    
    %% Main Output
    
    % How many stem cells at each age have a pathogenic mutation present.
    
    for uu = 1 : gg.numDiv
        Mut = find(MutatedAll(uu,:) >= (gg.mtDNA*gg.mutThreshold));
        MutNo = numel(Mut);
        MutatedSCAge(uu,pp) = MutNo;
    end
    
    %% To match the biological data I need to identify the clonally
    % expanded mutations (>25% heteroplasmy) at 70 years of age
    % equivalent to 3647 numDivs.
    
    for cc = 1  : gg.initS
        
        speciesPresent(:,cc) = mtDNAmutations(3647,((gg.mtDNA*cc)-(gg.mtDNA-1)):(gg.mtDNA*cc))';
        
    end
    
    % speciesPresent now gives the mutation
    % For each stem cell,find unique values and see if any of them are over 25%
    
    SingleMut = zeros(1,gg.initS);
    MultipleMut = zeros(1,gg.initS);
    
    for cc = 1 : gg.initS
        
        % Clonally expanded point mutation present?
        
        temp = unique(speciesPresent(:,cc));
        temp(temp==0) = [];
        
        % temp contains all the mtDNA mutations that are present at the age
        % of 70 years. Need to know if any of these are present in
        % speciesIDRecord.
        
        if isempty(temp)
            
        else
            
            for xx = 1 : numel(temp)
                
                temp2 = numel(find(speciesPresent(:,cc) == temp(xx)));
                temp3 = temp2 / gg.mtDNA * 100;
                
                findDouble = find(speciesIDRecord == temp(xx));
                
                if temp3 > 25 && isempty(findDouble)
                    SingleMut(cc) = SingleMut(cc) + 1;
                    MultipleMut(cc) = MultipleMut(cc) + 1;
                    
                elseif temp3 > 25 && ~isempty(findDouble)
                    MultipleMut(cc) = MultipleMut(cc) + speciesIDMultRecord(findDouble);
                    
                end
                
            end
            
        end
        
    end
    
    SingleMutRecord(pp,:) = SingleMut;
    MultipleMutRecord(pp,:) = MultipleMut;
    
    % Work out the probability for each age (numDivs) that there will be a
    % mutation present
    
    for tt = 1 : gg.initS
        for ii = 1 : gg.numDiv
            mutProbAge(2,ii) = mutProbAge(2,ii) + 1;
            if MutatedAll(ii,tt) >= gg.mtDNA*gg.mutThreshold
                mutProbAge(1,ii) = mutProbAge(1,ii) + 1;
            end
        end
    end
    
    % Update waitbar
    
    waitbar(pp/gg.numRuns,h)
    
end

delete(h)

gg.FailedCE(:,1) = [];
gg.SuccessCE(:,1) = [];

FailedClonal = gg.FailedCE;
SuccessClonal = gg.SuccessCE;


gg.NicheFailedSC(:,1) = [];
gg.NicheSuccessSC(:,1) = [];

NicheFailed = gg.NicheFailedSC;
NicheSuccess = gg.NicheSuccessSC;

% Make single and multiple mutation records spit out similar data to
% Biological Data

a = find(SingleMutRecord > 0);
a1 = find(MultipleMutRecord > 0);

b = numel(a);
b1 = numel(a1);

for kk = 1 : 20
    
    c = find(SingleMutRecord == kk);
    c1 = find(MultipleMutRecord == kk);
    
    SingleMutRecordResult(kk) = (numel(c)) / b * 100;
    MultipleMutRecordResult(kk) = (numel(c1)) / b1 * 100;
    
end

% Save the single and multiple mtDNA mutation data

gg.SingleMutRecordResult = SingleMutRecordResult;
gg.MultipleMutRecordResult = MultipleMutRecordResult;

% Save the mutation probabilities by age

mutProbAgeFinal = mutProbAge(1,:) ./ mutProbAge(2,:) * 100;

gg.mutProbAgeFinal = mutProbAgeFinal;

end