function zombiesInDenmark4()
%% Exact stochastic Simulation - Direct Method (Gillespie)

% 1. Initialize
inp = load('DenmarkMapWithInhab.mat');
inp2 = load('adjacencyMatrix.mat');
adjMat = inp2.spar_adjacency_matrix;
%numAdj = nnz(adjMat);
%[rowsAdj, colsAdj] = find(adjMat);

data = inp.data;
numInhab0 = [data{4,:}];
% define
% state = [S1, S2, S3,...
%          E1, E2, E3,...
%          Z1, Z2, Z3,...
%          R1, R2, R3,...]


numZones = length(numInhab0);

state = zeros(4,numZones);
state(1,:) = numInhab0;

% Define initial zombies
NZom = 10; % Number of zombies in initial state
zomIndex = 1; % The land area, that has the initial zombies.
state(3, zomIndex) = NZom;


t = 0;
tFinal = 100000;
n = 1;
nMax = 100000;

% Rates
beta = 3.6e-4;% /hr/person %Bite
alpha = 0.8; % inverse virulence
kappa = alpha*beta;% /hr/person % Kill
nu = 2;%/hr % Transform
mu = 0.914;% /hr Move
 

% REACTIONS
%a) Si + Zi -> Si + Ei
%b) Ei -> Zi
%c) Si + Zi -> Si + Ri
%d) Zi -> Zj for all adjacent areas

% Save state transitions 
n_transitions_per_file = 10000;

transitions = zeros(4,n_transitions_per_file);
transitionsFileName = 'zomInDKrun2.mat';
res.beta = beta;
res.alpha = alpha;
res.nu = nu;
res.mu = mu;
res.initState = state;
save(transitionsFileName, 'res')
% mapping of REACTION to int 
%) a) -> 1
%) b) -> 2
%) c) -> 3
%) d) -> 4



while t  < tFinal &&  n <= nMax
    %2. Calculate the propensity function for all i
    nonZeroFilter = state(1:3,:) > 0;
    
    SZFilter = nonZeroFilter(1,:) & nonZeroFilter(3,:);
    aa = beta*state(1,SZFilter).*state(3,SZFilter);
    len_aa = length(aa);
    ab = nu*state(2,nonZeroFilter(2,:));
    len_ab = length(ab);
    ac = alpha*aa;
    len_ac = len_aa;
    % Find adjacencies between zones that have zombies in them
    [adjWithZ, adjWithZ2] = find(adjMat & nonZeroFilter(3,:).'); % Zombies can potentially
    % move across any adjacency and each adjacency has two directions. However,
    %we only care about the directions where there are zombies to move.
    [adjWithZBack, adjWithZBack2] = find(adjMat & nonZeroFilter(3,:));
    ad = mu*[state(3,adjWithZ) state(3,adjWithZBack2)];
    
    if len_aa + len_ab + len_ac + length(ad) == 0
        disp('The zombies in Denmark have been defeated!')
        score = sum(state,2);
        fprintf(['Stats:\n\tSusceptibles: ' num2str(score(1)) '\n\tExposed: ' num2str(score(2)) '\n\tZombies: ' num2str(score(3)) '\n\tRemoved: ' num2str(score(4)) '\n' ])
        
        break
    end
    
    a = [aa ab ac ad];
    %3. For each i, generate a putative time, ti, according to an exponential distribution with parameter ai. 
    ts = exprnd(1./a);
    %pause()
    
    %4. Let i_min be the reaction for which ti is least;
    [t_min, i_min] = min(ts);
    
    % NOTE: 
    %   i_min is an index ts
    %   shifting i_min by len_XX makes makes it an index in
    %   state(1,nonZero), call it i_min'
    %   finding the index of i_min' in nonZeroFilter, gives the index of
    %   the zone where the change happens, ind
    
    %6. Change the number of molecules to reflect execution of reaction.
    if i_min <= len_aa % perform bite
        tmp = find(SZFilter,i_min,'first');
        ind = tmp(end);
        %disp(['Human bit in ' num2str(i_min)])
        state(1:2,ind) = state(1:2,ind) + [-1; 1];
        trans = [1, ind, 0];
    elseif i_min <= len_aa + len_ab % perform transform
        i_minp = i_min - len_aa;
        tmp = find(nonZeroFilter(2,:),i_minp,'first');
        ind = tmp(end);
        %disp(['Zombie transformed in ' num2str(ind)])
        state(2:3,ind) = state(2:3,ind) + [-1; 1];
        trans = [2, ind, 0];
    elseif i_min <= len_aa + len_ab + len_ac % perform kill
        i_minp = i_min - (len_aa + len_ab);
        tmp = find(SZFilter,i_minp,'first');
        ind = tmp(end);
        %disp(['Zombie killed in ' num2str(ind)])
        state(3:4,ind) = state(3:4,ind) + [-1; 1];
        trans = [3, ind, 0];
    elseif i_min <= len_aa + len_ab + len_ac + length(adjWithZ) % perform move A -> B
        ind = i_min - (len_aa + len_ab + len_ac); % the index in "rows"/"cols"  corresponding to the move
        %disp(['Zombie moved from ' num2str(adjWithZ(ind)) ' to ' num2str(adjWithZ2(ind)) ])
        state(3,[adjWithZ(ind) adjWithZ2(ind)]) = state(3,[adjWithZ(ind) adjWithZ2(ind)]) + [-1 1];
        trans = [4, adjWithZ(ind), adjWithZ2(ind)];
    else % perform move B -> A
        ind = i_min - (len_aa + len_ab + len_ac + length(adjWithZ)); % the index in "rows"/"cols"  corresponding to the move
        %disp(['Zombie moved from ' num2str(adjWithZBack2(ind)) ' to ' num2str(adjWithZBack(ind)) ])
        state(3,[adjWithZBack2(ind) adjWithZBack(ind)]) = state(3,[adjWithZBack2(ind) adjWithZBack(ind)]) + [-1 1];
        trans = [4, adjWithZBack2(ind), adjWithZBack(ind)];
    end
    np = mod(n,n_transitions_per_file);
    t = t + t_min;
    
    if  np ~= 0
        transitions(:,np) = [t trans].';
    else
        transitions(:,n_transitions_per_file) = [t trans].';
        res.(['transitions' num2str(n/n_transitions_per_file)]) = transitions;
        res.finalState = state;
        save(transitionsFileName, 'res', '-append')
        disp('')
        %inp = load(transitionsFileName, 'res');
        %pause()
    end
    n = n + 1;
end

