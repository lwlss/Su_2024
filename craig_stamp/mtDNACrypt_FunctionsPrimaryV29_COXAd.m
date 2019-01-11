%% Niche Succession Model - FINAL
% Stem cell dynamics and mutated mtDNA clonal expansion
%
% The script is an amalgamation of the previous crypt model. It identifies
% that there are a certain number of mtDNA molecules residing within each
% stem cell of the crypt. With the evolution of stem cell divisions, the
% number of mutated mtDNA molecules evolves stochastically according to
% pre-determined probabilities. Also, with each additional mutated mtDNA
% molecule, the model determines which kind of mutation has developed
% according to probability data previously aquired. Therefore, this model
% is a more accurate representation of the processes that take place within
% the crypt and at the tissue level.

% v8 v7.3 compression of the saved variables. Time bar for each crypt
% simulation
% v9 Intergrates a user interface box which asks all of the required
% parameters
% v10 enables the user to open and close a parameter list file while the
% simulation is running to add new simulations to the list
% v11 Crypt fission fix which allows the contiuation of the script for
% large numbers of runs
% v12 The way in which the new crypt data from fission is intergrated into
% the final results table is drawn using the randperm function therefore
% unique numbers are now selected.
% v13 Every new crypt fission event does not overwrite the resultant data
% from a previous crypt result with the addition of a cryptReplaceCount
% Counter for every crypt that is replace in the final data set.
% v14 Solving memory issues which arise after prolonged model simulation,
% solved by executing clear command before everyne model of a different
% parameter set.
% v26 Fully COX deficienct stem cells are subject to random removal for some
% relevent biological reason -- at the mtDNA molecular level. SpeciesID
% identification and removal.
% v27 Additional COX deficient stem cell division and ParameterNames
% variable has been transposed.
% v28 Mitochondrial degradation incorporated into the model
% v29 Intergration of the new transition matrices that take into account
% the possible asymmetric segregation of mutated mtDNA molecules.

% Load the parameter list file
load ParameterListFittingScan

% Determine how many simulations are to be carried out using 'cycle'

ParameterNames = ParameterNames';
a = size(ParameterNames);
b = a(1);
cycle = b - 1;

% Set 'qq' to 1 for the first simulation
qq = 1;

% For every simulation increase 'qq' by 1 until total number of simulations
% 'cycle' has been reached

while qq <= cycle
    
    % Clear memory after every simulation so that the memory doesnt become too
    % fragmented when many simulations are to be carried out
    save SystemMemoryClearUp qq
    save SystemMemoryClearUp cycle -append
    clear
    
    % Set global varable structure 'gg' where all parameters are stored for the
    % model simulation and where all metrics are stored once model is completed
    global gg
    global dd
    
    % Reload all critical variables after the memory purge
    load SystemMemoryClearUp
    load ParameterListFittingScan
    
    % Transpose parameter variables so that they're in the correct format
    ParameterNames = ParameterNames';
    
    % What is the filename for the overall results?
    gg.finalFilename = datestr(clock,30);
    ParameterNames(qq+1, 14) = {gg.finalFilename};
    
    % Shuffle random number generator before every simulation so that the model
    % is truly stochastic in nature
    rng('shuffle');
    
    %% Load all the variables into the 'gg' global variable
    
    % Number of crypts generated per simulation
    gg.numRuns = cell2mat(ParameterNames(qq+1,1));
    
    % The percentage threshold that characterises a stem cell as COX deficient
    gg.mutThreshold = cell2mat(ParameterNames(qq+1,2));
    
    % The number of asynchronous stem cell divisions that portrays the human
    % lifespan (1 stem cell division per week)
    gg.numDiv = cell2mat(ParameterNames(qq+1,3));
    
    % Number of mtDNA molecules contained within each stem cell
    gg.mtDNA = cell2mat(ParameterNames(qq+1,4));
    
    % Number of stem cells contained within crypts
    gg.initS = cell2mat(ParameterNames(qq+1,5));
    
    % Stem cell division types 'Pa' Asymmetric probability 'Ps' Symmetric
    % probability
    gg.Pa = cell2mat(ParameterNames(qq+1,6));
    gg.Ps = (1 - gg.Pa)/2;
    
    % Advantage to COX deficient stem cell to divide more often according to
    % 'adv' which increases Ps and reduces Pa1
    gg.adv = cell2mat(ParameterNames(qq+1,7));
    
    % Which method will be used for the mutation rate?
    gg.mutMethod = char(ParameterNames(qq+1,8));
    
    % What is the base mutation rate?
    gg.mutationRate = cell2mat(ParameterNames(qq+1,9));
    
    % Is there a mutation rate fold change from 0 to 80 years of age? if so
    % what is it? If there isn't this should be set to 1.
    gg.mutationRateFold = cell2mat(ParameterNames(qq+1,15));
    
    % Calculating the mutation rate vector for each division step in the model
    c = gg.mutationRate;
    m = (gg.mutationRateFold *gg.mutationRate - gg.mutationRate) / 4171;
    for zz = 1 : 5211
        gg.mutationRate1(zz) = m*zz + c;
    end
    
    %% COX Correction Factors
    
    gg.COXCorrectionFactor = cell2mat(ParameterNames(qq+1,16));
    
    gg.COXCorrectionFactor2 = cell2mat(ParameterNames(qq+1,17));
    
    gg.COXSCTimePoint = char(ParameterNames(qq+1,18));
    
    gg.COXSCTimePointInterval = cell2mat(ParameterNames(qq+1,19));
    
    gg.COXDefCycleRepeats = cell2mat(ParameterNames(qq+1,20));
    
    gg.MitoDegradation = cell2mat(ParameterNames(qq+1,21));
    
    %% Crypt Fission
    
    gg.cryptFission = char(ParameterNames(qq+1,11));
    
    if strcmp('yes',gg.cryptFission)
        
        gg.cryptNormalPercentage = cell2mat(ParameterNames(qq+1,12));
        gg.cryptFissionFactor = cell2mat(ParameterNames(qq+1,13));
        gg.cryptFissionProb = (1/gg.numDiv)*gg.cryptNormalPercentage;
        gg.cryptFisSave = 1;
        
    end
    
    % Parameters to prime the metrics to be recorded
    
    gg.FailedCE = [0;0;0];
    gg.SuccessCE = [0;0;0];
    
    gg.NicheFailedSC = [0;0;0];
    gg.NicheSuccessSC = [0;0;0];
    
    % Load probability tables and mutation probabilities
    
    switch gg.mtDNA
        
        case 5
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb5.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb5.mat');
            
        case 10
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb10.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb10.mat');
            
        case 25
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb25.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb25.mat');
            
        case 50
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb50.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb50.mat');
            
        case 100
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb100.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100.mat');
            
            % Load the advantage DivProbs for 100 mtDNA SCs
            
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb10010.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100100.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb1002.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb10011.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100101.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100102.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100103.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb100108.mat');
            % load('D:\Niche Succession Model Transfer\dividingMutations\DivProb1001001.mat');
            
        case 200
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb200.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb200.mat');
            
        case 400
            
            load('D:\Niche Succession Model Transfer\replicatingMutations\RepProb400.mat');
            load('D:\Niche Succession Model Transfer\dividingMutations\DivProb400.mat');
            
        otherwise
            warning('Please enter a valid mtDNA number for which a transition matrix has been created.');
            
    end
    
    gg.RepProb = RepProb; clearvars RepProb
    gg.DivProb = DivProb; clearvars DivProb
    
    % gg.DivProb10 = DivProb10010; clearvars DivProb10010
    % gg.DivProb100 = DivProb100100; clearvars DivProb100100
    % gg.DivProb2 = DivProb1002; clearvars DivProb1002
    % gg.DivProb11 = DivProb10011; clearvars DivProb10011
    
    % gg.DivProb101 = DivProb100101; clearvars DivProb100101
    % gg.DivProb102 = DivProb100102; clearvars DivProb100102
    gg.DivProb103 = DivProb100103; clearvars DivProb100103
    % gg.DivProb108 = DivProb100108; clearvars DivProb100108
    % gg.DivProb1001 = DivProb100108; clearvars DivProb100108
    
    % Least sqaures determination
    
    gg.LeastSqauresRunInterval = gg.numRuns / 100;
    
    % Save parameters name with new information
    
    ParameterNames = ParameterNames';
    
    save ParameterListFittingScan ParameterNames
    
    clearvars ParameterNames
    
    tic
    
    if strcmp('yes',gg.cryptFission)
        
        % Which simulation type is going to be performed
        
        switch gg.mutMethod
            
            case 'constant'
                
                [MutatedSCAgeFinal, MutatedSCAgeCorrFinal,...
                    MutatedSCAgeCorr2Final] = mtDNACrypt_ConstantV11FCN_COXAd_CF();
                
            case 'exponential'
                
                MutatedSCAgeFinal = mtDNACrypt_ExponentialV2FCN_CF();
        end
        
        save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
        save MutatedSCAgeCorrected MutatedSCAgeCorrFinal -v7.3
        save MutatedSCAgeCorrected2 MutatedSCAgeCorr2Final -v7.3
        
        % Load the crypt data where crypt fission has occurred
        
        rr = load('cryptFissionResult.mat');
        rr = rmfield(rr,'Kickstart1');
        cryptNames = fieldnames(rr);
        
        % The loaded crypt data is loaded into a cell in order to find out the
        % number of crypts that underwent fission
        
        s = numel(cryptNames);
        
        % Record how many doublets have formed
        
        gg.colonys.doublets = s;
        
        % Bring out all the crypts data before crypt fission occurred from the
        % structured array
        
        struct2var(rr);
        
        % Create a new matrix that contains the continued data from the crypt
        % fission event
        
        MutatedSCAgeFission = zeros(gg.numDiv,s);
        
        % Create a new vector which records the point at which each crypt fission
        % event took place so this can be used only import new crypt fission data
        % into the original results matrix
        
        CryptFisTime = zeros(1,s);
        
        % Reset the filename so that the new crypt fission data is saved to a
        % different location
        
        gg.filename = 'cryptFissionResult2';
        
        h = waitbar(0,'Simulating model with crypt fission - part 2, please wait...');
        
        for kk = 1 : s
            
            s1 = cryptNames(kk,1);   % Identify a single crypt fission event
            s2 = char(s1);  % Convert the crypt name into a string
            s3 = eval(s2);  % Assign the matrix to the string
            
            
            if strcmp('constant',gg.mutMethod)
                [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                    mtDNACrypt_ConstantV2FCN_CF_Input(s3);
            else
                [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                    mtDNACrypt_ExponentialV2FCN_CF_Input(s3);
            end
            
            waitbar(kk/s,h)
            
        end
        
        delete(h)
        
        % Save the crypts that have undergone a second round of fission
        
        save(gg.filename,'dd','-v7.3');
        
        gg.FissionTime = CryptFisTime;
        
        load MutatedSCAgeFission MutatedSCAgeFinal
        
        % Replace random crypts
        
        cryptReplaceCount = 0;
        
        a = numel(CryptFisTime);
        b = size(MutatedSCAgeFinal);
        cryptReplaceNo = [randperm(b(2)) randperm(b(2)) randperm(b(2))];
        
        %     c = ceil(rand(1,a)*b(2));
        
        for tt = 1 : a
            MutatedSCAgeFinal(CryptFisTime(tt):gg.numDiv,cryptReplaceNo(tt + cryptReplaceCount)) =...
                MutatedSCAgeFission(CryptFisTime(tt):gg.numDiv,tt);
        end
        
        cryptReplaceCount = cryptReplaceCount + a;
        
        save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
        
        clearvars -except gg cycle qq cryptReplaceCount cryptReplaceNo
        
        gg.analysisComplete = 0;
        
        %% Analysis OR Crypt Fission Round 2
        
        load('cryptFissionResult2.mat');
        rr = dd;
        clearvars dd
        
        if isstruct(rr)
            cryptNames = fieldnames(rr);
        else
            cryptNames = [];
        end
        
        if isempty(cryptNames)
            moreCryptFission = 'no';
        else
            moreCryptFission = 'yes';
        end
        
        if strcmp('no',moreCryptFission)
            
            gg.analysisComplete = 1;
            
        else
            
            % The loaded crypt data is loaded into a cell in order to find out the
            % number of crypts that underwent fission
            
            s = numel(cryptNames);
            
            % Record how many doublets have formed
            
            gg.colonys.triplets = s;
            
            % Bring out all the crypts data before crypt fission occurred from the
            % structured array
            
            struct2var(rr);
            
            % Create a new matrix that contains the continued data from the crypt
            % fission event
            
            MutatedSCAgeFission = zeros(gg.numDiv,s);
            
            % Create a new vector which records the point at which each crypt fission
            % event took place so this can be used only import new crypt fission data
            % into the original results matrix
            
            CryptFisTime = zeros(1,s);
            
            % Reset the filename so that the new crypt fission data is saved to a
            % different location
            
            gg.filename = 'cryptFissionResult3';
            
            clearvars -global dd
            global dd
            
            h = waitbar(0,'Simulating model with crypt fission - part 3, please wait...');
            
            for kk = 1 : s
                
                s1 = cryptNames(kk,1);   % Identify a single crypt fission event
                s2 = char(s1);  % Convert the crypt name into a string
                s3 = eval(s2);  % Assign the matrix to the string
                
                
                if strcmp('constant',gg.mutMethod)
                    [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                        mtDNACrypt_ConstantV2FCN_CF_Input(s3);
                else
                    [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                        mtDNACrypt_ExponentialV2FCN_CF_Input(s3);
                end
                
                waitbar(kk/s,h)
                
            end
            
            delete(h)
            
            save(gg.filename,'dd','-v7.3')
            
            gg.FissionTime = [gg.FissionTime CryptFisTime];
            
            load MutatedSCAgeFission MutatedSCAgeFinal
            
            % Replace random crypts
            
            a = numel(CryptFisTime);
            
            for tt = 1 : a
                MutatedSCAgeFinal(CryptFisTime(tt):gg.numDiv,cryptReplaceNo(tt + cryptReplaceCount)) =...
                    MutatedSCAgeFission(CryptFisTime(tt):gg.numDiv,tt);
            end
            
            cryptReplaceCount = cryptReplaceCount + a;
            
            save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
            
            clearvars -except gg cycle qq cryptReplaceCount cryptReplaceNo
            
        end
        
        %% Analysis OR Crypt Fission Round 3
        
        if gg.analysisComplete ~= 1;
            
            load('cryptFissionResult3.mat');
            rr = dd;
            clearvars dd
            
            if isstruct(rr)
                cryptNames = fieldnames(rr);
            else
                cryptNames = [];
            end
            
            if isempty(cryptNames)
                moreCryptFission = 'no';
            else
                moreCryptFission = 'yes';
            end
            
            if strcmp('no',moreCryptFission)
                
                gg.analysisComplete = 1;
                
            else
                
                % The loaded crypt data is loaded into a cell in order to find out
                % the number of crypts that underwent fission
                
                s = numel(cryptNames);
                
                % Bring out all the crypts data before crypt fission occurred from
                % the structured array
                
                struct2var(rr);
                
                % Record how many doublets have formed
                
                gg.colonys.quadruplets = s;
                
                % Create a new matrix that contains the continued data from the
                % crypt fission event
                
                MutatedSCAgeFission = zeros(gg.numDiv,s);
                
                % Create a new vector which records the point at which each crypt
                % fission event took place so this can be used only import new
                % crypt fission data into the original results matrix
                
                CryptFisTime = zeros(1,s);
                
                % Reset the filename so that the new crypt fission data is savd to
                % a different location
                
                gg.filename = 'cryptFissionResult4';
                
                clearvars -global dd
                global dd
                
                h = waitbar(0,'Simulating model with crypt fission - part 4, please wait...');
                
                for kk = 1 : s
                    
                    s1 = cryptNames(kk,1); % Identify a single crypt fission event
                    s2 = char(s1);  % Convert the crypt name into a string
                    s3 = eval(s2);  % Assign the matrix to the string
                    
                    if strcmp('constant',gg.mutMethod)
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ConstantV2FCN_CF_Input(s3);
                    else
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ExponentialV2FCN_CF_Input(s3);
                    end
                    
                    waitbar(kk/s,h)
                    
                end
                
                delete(h)
                
                save(gg.filename,'dd','-v7.3');
                
                gg.FissionTime = [gg.FissionTime CryptFisTime];
                
                load MutatedSCAgeFission MutatedSCAgeFinal
                
                % Replace random crypts
                
                a = numel(CryptFisTime);
                
                for tt = 1 : a
                    MutatedSCAgeFinal(CryptFisTime(tt):5210,cryptReplaceNo(tt + cryptReplaceCount)) =...
                        MutatedSCAgeFission(CryptFisTime(tt):5210,tt);
                end
                
                cryptReplaceCount = cryptReplaceCount + a;
                
                save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
                
                clearvars -except gg cycle qq cryptReplaceCount cryptReplaceNo
                
            end
            
        end
        
        %% Analysis OR Crypt Fission Round 4
        
        if gg.analysisComplete ~= 1;
            
            load('cryptFissionResult4.mat');
            rr = dd;
            clearvars dd
            
            if isstruct(rr)
                cryptNames = fieldnames(rr);
            else
                cryptNames = [];
            end
            
            if isempty(cryptNames)
                moreCryptFission = 'no';
            else
                moreCryptFission = 'yes';
            end
            
            if strcmp('no',moreCryptFission)
                
                gg.analysisComplete = 1;
                
            else
                
                % The loaded crypt data is loaded into a cell in order to find out
                % the number of crypts that underwent fission
                
                s = numel(cryptNames);
                
                % Record how many doublets have formed
                
                gg.colonys.quintuplets = s;
                
                % Bring out all the crypts data before crypt fission occurred from
                % the structured array
                
                struct2var(rr);
                
                % Create a new matrix that contains the continued data from the
                % crypt fission event
                
                MutatedSCAgeFission = zeros(gg.numDiv,s);
                
                % Create a new vector which records the point at which each crypt
                % fission event took place so this can be used only import new
                % crypt fission data into the original results matrix
                
                CryptFisTime = zeros(1,s);
                
                % Reset the filename so that the new crypt fission data is savd to
                % a different location
                
                gg.filename = 'cryptFissionResult5';
                
                clearvars -global dd
                global dd
                
                h = waitbar(0,'Simulating model with crypt fission - part 5, please wait...');
                
                for kk = 1 : s
                    
                    s1 = cryptNames(kk,1);  % Identify a single crypt fission event
                    s2 = char(s1);  % Convert the crypt name into a string
                    s3 = eval(s2);  % Assign the matrix to the string
                    
                    if strcmp('constant',gg.mutMethod)
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ConstantV2FCN_CF_Input(s3);
                    else
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ExponentialV2FCN_CF_Input(s3);
                    end
                    
                    waitbar(kk/s,h)
                    
                end
                
                delete(h)
                
                save(gg.filename,'dd','-v7.3');
                
                gg.FissionTime = [gg.FissionTime CryptFisTime];
                
                load MutatedSCAgeFission MutatedSCAgeFinal
                
                % Replace random crypts
                
                a = numel(CryptFisTime);
                %             b = size(MutatedSCAgeFinal);
                %             c = randperm(b(2),a);
                
                %             c = ceil(rand(1,a)*b(2));
                
                for tt = 1 : a
                    MutatedSCAgeFinal(CryptFisTime(tt):5210,cryptReplaceNo(tt + cryptReplaceCount)) =...
                        MutatedSCAgeFission(CryptFisTime(tt):5210,tt);
                end
                
                cryptReplaceCount = cryptReplaceCount + a;
                
                save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
                
                clearvars -except cycle qq gg cryptReplaceCount cryptReplaceNo
                
            end
            
        end
        
        %% Analysis OR Crypt Fission Round 5 (Last Round)
        
        if gg.analysisComplete ~= 1;
            
            load('cryptFissionResult5.mat');
            rr = dd;
            clearvars dd
            
            if isstruct(rr)
                cryptNames = fieldnames(rr);
            else
                cryptNames = [];
            end
            
            if isempty(cryptNames)
                moreCryptFission = 'no';
            else
                moreCryptFission = 'yes';
            end
            
            if strcmp('no',moreCryptFission)
                
                gg.analysisComplete = 1;
                
            else
                
                % The loaded crypt data is loaded into a cell in order to find out
                % the number of crypts that underwent fission
                
                s = numel(cryptNames);
                
                % Record how many doublets have formed
                
                gg.colonys.sextuplets = s;
                
                % Bring out all the crypts data before crypt fission occurred from
                % the structured array
                
                struct2var(rr);
                
                % Create a new matrix that contains the continued data from the
                % crypt fission event
                
                MutatedSCAgeFission = zeros(gg.numDiv,s);
                
                % Create a new vector which records the point at which each crypt
                % fission event took place so this can be used only import new
                % crypt fission data into the original results matrix
                
                CryptFisTime = zeros(1,s);
                
                % Reset the filename so that the new crypt fission data is savd to
                % a different location
                
                gg.filename = 'cryptFissionResult6';
                
                clearvars -global dd
                global dd
                
                h = waitbar(0,'Simulating model with crypt fission - part 6, please wait...');
                
                for kk = 1 : s
                    
                    s1 = cryptNames(kk,1);  % Identify a single crypt fission event
                    s2 = char(s1);  % Convert the crypt name into a string
                    s3 = eval(s2);  % Assign the matrix to the string
                    
                    if strcmp('constant',gg.mutMethod)
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ConstantV2FCN_CF_Input(s3);
                    else
                        [MutatedSCAgeFission(:,kk),CryptFisTime(kk)] =...
                            mtDNACrypt_ExponentialV2FCN_CF_Input(s3);
                    end
                    
                    waitbar(kk/s,h)
                    
                end
                
                delete(h)
                
                save(gg.filename,'dd','-v7.3');
                
                gg.FissionTime = [gg.FissionTime CryptFisTime];
                
                load MutatedSCAgeFission MutatedSCAgeFinal
                
                % Replace random crypts
                
                a = numel(CryptFisTime);
                
                for tt = 1 : a
                    MutatedSCAgeFinal(CryptFisTime(tt):5210,cryptReplaceNo(tt + cryptReplaceCount)) =...
                        MutatedSCAgeFission(CryptFisTime(tt):5210,tt);
                end
                
                cryptReplaceCount = cryptReplaceCount + a;
                
                save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
                
                clearvars -except cycle qq gg cryptReplaceCount cryptReplaceNo
                
            end
            
        end
        
        %% What happens if there is no crypt fission that occurs??
        
    else
        
        switch gg.mutMethod
            
            case 'constant'
                
                [MutatedSCAgeFinal, MutatedSCAgeCorrFinal, MutatedSCAgeCorr2Final] = mtDNACrypt_ConstantV11FCN_COXAd();
                
            case 'exponential'
                
                MutatedSCAgeFinal = mtDNACrypt_ExponentialV2FCN();
        end
        
        save MutatedSCAgeFission MutatedSCAgeFinal -v7.3
        save MutatedSCAgeCorrected MutatedSCAgeCorrFinal -v7.3
        save MutatedSCAgeCorrected2 MutatedSCAgeCorr2Final -v7.3
        
        
        clearvars -except gg cycle qq
        
    end
    
    %% Least squares to determine the optimum number of runs
    
    load MutatedSCAgeFission MutatedSCAgeFinal
    load MutatedSCAgeCorrected MutatedSCAgeCorrFinal
    load MutatedSCAgeCorrected2 MutatedSCAgeCorr2Final
    
    for jj = 1 : gg.numRuns / gg.LeastSqauresRunInterval
        for ii = 1 : gg.numDiv
            a(ii,jj) = mean(MutatedSCAgeFinal(ii,1:(jj*gg.LeastSqauresRunInterval)));
        end
    end
    
    b = size(a);
    
    for tt = 1 : b(2)
        c(tt) = sum(a(:,tt));
    end
    
    for rr = 1 : length(c)-1
        d(rr) = sqrt((c(rr+1) - c(rr))^2);
    end
    
    
    %% Save Simulation Results
    
    % remove fields that are no longer required
    
    if strcmp('no',gg.cryptFission)
        fields = {'RepProb','DivProb'};
    else
        fields = {'filename','analysisComplete','RepProb','DivProb','cryptFisSave'};
    end
    
    gg = rmfield(gg,fields);
    
    % Make the final results and the least square values global gg
    
    gg.MutatedSCAgeFinal = MutatedSCAgeFinal;
    gg.MutatedSCAgeFinalCorr = MutatedSCAgeCorrFinal;
    gg.MutatedSCAgeFinalCorr2 = MutatedSCAgeCorr2Final;
    
    gg.LeastSqaures = d;
    
    % save the whole global gg variable
    
    gg.SimulationTime = toc;
    
    save(gg.finalFilename,'gg','-v7.3')
    
    % Delete files that are no longer required
    
    if strcmp('no',gg.cryptFission)
        delete('MutatedSCAgeFission.mat')
        delete('MutatedSCAgeCorrected.mat')
        delete('MutatedSCAgeCorrected2.mat')
        delete('SystemMemoryClearUp.mat')
    else
        delete('MutatedSCAgeFission.mat')
        delete('cryptFissionResult*')
        delete('MutatedSCAgeCorrected.mat')
        delete('MutatedSCAgeCorrected2.mat')
        delete('SystemMemoryClearUp.mat')
    end
    
    % Clear up the desktop and workspace and declare the model is finished
    
    qqstr = num2str(qq);
    disp(['Model ' qqstr ' Completed Successfully']);
    
    % Assess the number of parameters that
    
    load ParameterListFittingScan
    ParameterNames = ParameterNames';
    a = size(ParameterNames);
    cycle = a(1) - 1;
    clearvars -except qq cycle
    
    qq = qq + 1;
    
end

disp('Simulations Complete');
h = msgbox('Simulation Complete','Success');
clear all;

    


