function zombiesInDenmark4VOLSOM()
%% Exact stochastic Simulation - Direct Method (Gillespie)

% 1. Initialize
inp = load('data/DenmarkMapWithInhab.mat');
inp2 = load('data/adjacencyMatrix.mat');
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

%inp = load('../data/zomInDKrun1FinalState.mat');
state = zeros(5,numZones);
state(1,:) = numInhab0;


%state = zeros(4,numZones);
%state(1,:) = numInhab0;

% Define initial zombies
NZom = 10; % Number of zombies in initial state
zomIndex = 1357; %%% Skagen%%% The land area, that has the initial zombies.
state(3, zomIndex) = NZom;

% Define initial military
NMilitary = 2.1e4; 
inp = load('data/baseParishes.mat');
% Indices of parishes that have military bases
base_indices = inp.base_parishes;
% Number of these
num_bases = length(base_indices);
mili_per_base = floor(NMilitary/num_bases);
% Assign initial personel
state(5,base_indices) = mili_per_base;


t = 0;
n = 1;
nMax = 8e7;
% Rates
beta = 2.87e-4;% /hr/person %Bite
alpha = 0.8; % inverse virulence
kappa = alpha*beta;% /hr/person % Kill
nu = 2;%/hr % Transform
mu = 0.0568;% /hr Move
tau = 0.0023;% /hr Train 
f = 10;
mu_M = f*mu;% /hr Move military
 

% REACTIONS
%a) Si + Zi -> Zi + Ei
%b) Ei -> Zi
%c) Si + Zi -> Si + Ri
%d) Zi -> Zj for all adjacent areas
%e) Si -> Mi 
%f) Mi + Zi -> Mi + Ri
%g) Mi + Zi -> Ei + Zi 
%h) Mi -> Mj for all adjacent areas


% Save state transitions 
n_transitions_per_file = 50000;

transitions = zeros(4,n_transitions_per_file);
transitionsFileName = '../data/zomInDKrun3.mat';
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
%) e) -> 5
%) f) -> 6
%) g) -> 7
%) h) -> 8

all_lengths = zeros(1, 10);

%try
while n <= nMax
    %2. Calculate the propensity function for all i
    nonZeroFilter = state > 0;
    
    
    SZFilter = nonZeroFilter(1,:) & nonZeroFilter(3,:);
    % a) Bite
    aa = beta*state(1,SZFilter).*state(3,SZFilter);
    all_lengths(1) = length(aa);
    % b) Transform
    ab = nu*state(2,nonZeroFilter(2,:));
    all_lengths(2) = length(ab);
    % c) Kill
    ac = alpha*aa;
    all_lengths(3) = all_lengths(1);
    % d) Move zom
    % Find adjacencies between zones that have zombies in them
    is_adj_has_ZMat = bsxfun(@and,adjMat,nonZeroFilter(3,:).');
    [adjWithZ, adjWithZ2] = find(is_adj_has_ZMat); % Zombies can potentially
    % move across any adjacency and each adjacency has two directions. However,
    %we only care about the directions where there are zombies to move.
    is_adj_has_ZMatBack = bsxfun(@and, adjMat, nonZeroFilter(3,:));
    [adjWithZBack, adjWithZBack2] = find(is_adj_has_ZMatBack);
    ad = mu*[state(3,adjWithZ) state(3,adjWithZBack2)];
    all_lengths(7:8) = [length(adjWithZ) length(adjWithZBack2)];
    % e) Train
    ae = tau*state(5,nonZeroFilter(5,:));
    all_lengths(4) = length(ae);
    % f) Bite M
    MZFilter = nonZeroFilter(5,:) & nonZeroFilter(3,:);
    af = beta/f*state(5,MZFilter).*state(3,MZFilter);
    all_lengths(5) = length(af);
    % g) Kill M
    ag = alpha*beta/f*state(5,MZFilter).*state(3,MZFilter);
    all_lengths(6) = all_lengths(5);
    % h) Move M
    [adjWithZAndM, adjWithZAndM2] = find(bsxfun(@and, is_adj_has_ZMat,nonZeroFilter(5,:).'));
    [adjWithZAndMBack, adjWithZAndMBack2] = find(bsxfun(@and, is_adj_has_ZMat,nonZeroFilter(5,:).'));
    ah = f*mu*[state(5,adjWithZAndM) state(5, adjWithZAndMBack2)];
    all_lengths(9:10) = [length(adjWithZAndM) length(adjWithZAndMBack2)];
    % Cumulated lengths
    cum_lengths = cumsum(all_lengths);
    if cum_lengths(end) == 0
        disp('The zombies in Denmark have been defeated!')
        score = sum(state,2);
        fprintf(['Stats:\n\tSusceptibles: ' num2str(score(1)) '\n\tExposed: ' num2str(score(2)) '\n\tZombies: ' num2str(score(3)) '\n\tRemoved: ' num2str(score(4)) '\n' ])
        
        break
    end
    
    %3. For each i, generate a putative time, ti, according to an exponential distribution with parameter ai. 
    ts = exprnd(1./[aa ab ac ae af ag ad ah]);
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
    if i_min <= cum_lengths(1) % perform bite
        tmp = find(SZFilter,i_min,'first');
        ind = tmp(end);
        %disp(['Human bit in ' num2str(i_min)])
        state(1:2,ind) = state(1:2,ind) + [-1; 1];
        trans = [1, ind, 0];
    elseif i_min <= cum_lengths(2) % perform transform
        i_minp = i_min - cum_lengths(1);
        tmp = find(nonZeroFilter(2,:),i_minp,'first');
        ind = tmp(end);
        %disp(['Zombie transformed in ' num2str(ind)])
        state(2:3,ind) = state(2:3,ind) + [-1; 1];
        trans = [2, ind, 0];
    elseif i_min <= cum_lengths(3) % perform kill
        i_minp = i_min - cum_lengths(2);
        tmp = find(SZFilter,i_minp,'first');
        ind = tmp(end);
        %disp(['Zombie killed in ' num2str(ind)])
        state(3:4,ind) = state(3:4,ind) + [-1; 1];
        trans = [3, ind, 0];
    elseif i_min <= cum_lengths(4) % Train military
        i_minp = i_min - cum_lengths(3);
        tmp = find(nonZeroFilter(5,:),i_minp,'first');
        ind = tmp(end);
        state([1 5],ind) = state([1 5],ind) + [-1; 1];
        trans = [5, ind, 0];
    elseif i_min <= cum_lengths(5)% Bite military
        i_minp = i_min - cum_lengths(4);
        tmp = find(MZFilter,i_minp,'first');
        ind = tmp(end);
        state([5 2],ind) = state([5 2],ind) + [-1; 1];
        trans = [6, ind, 0];
    elseif i_min <= cum_lengths(6)% Kill by military
        i_minp = i_min - cum_lengths(5);
        tmp = find(MZFilter,i_minp,'first');
        ind = tmp(end);
        %disp(['Zombie killed in ' num2str(ind)])
        state(3:4,ind) = state(3:4,ind) + [-1; 1];
        trans = [7, ind, 0];
    elseif i_min <= cum_lengths(7) % perform move A -> B
        ind = i_min - cum_lengths(6); % the index in "rows"/"cols"  corresponding to the move
        %disp(['Zombie moved from ' num2str(adjWithZ(ind)) ' to ' num2str(adjWithZ2(ind)) ])
        state(3,[adjWithZ(ind) adjWithZ2(ind)]) = state(3,[adjWithZ(ind) adjWithZ2(ind)]) + [-1 1];
        trans = [4, adjWithZ(ind), adjWithZ2(ind)];
    elseif i_min <= cum_lengths(8) % perform move B -> A
        ind = i_min - cum_lengths(7); % the index in "rows"/"cols"  corresponding to the move
        %disp(['Zombie moved from ' num2str(adjWithZBack2(ind)) ' to ' num2str(adjWithZBack(ind)) ])
        state(3,[adjWithZBack2(ind) adjWithZBack(ind)]) = state(3,[adjWithZBack2(ind) adjWithZBack(ind)]) + [-1 1];
        trans = [4, adjWithZBack2(ind), adjWithZBack(ind)];
    elseif i_min <= cum_lengths(9) % perform Military move A-> B
        ind = i_min - cum_lengths(8); % the index in "rows"/"cols"  corresponding to the move
        state(5,[adjWithZAndM(ind) adjWithZAndM2(ind)]) = state(5,[adjWithZAndM(ind) adjWithZAndM2(ind)]) + [-1 1];
        trans = [8, adjWithZAndM(ind), adjWithZAndM2(ind)];
    else % perform Military move B -> A
        ind = i_min - cum_lengths(9); % the index in "rows"/"cols"  corresponding to the move
        %disp(['Zombie moved from ' num2str(adjWithZBack2(ind)) ' to ' num2str(adjWithZBack(ind)) ])
        state(5,[adjWithZAndMBack2(ind) adjWithZAndMBack(ind)]) = state(5,[adjWithZAndMBack2(ind) adjWithZAndMBack(ind)]) + [-1 1];
        trans = [8, adjWithZAndMBack2(ind), adjWithZAndMBack(ind)];
    end
    np = mod(n,n_transitions_per_file);
    t = t + t_min;
    
    if  np ~= 0
        transitions(:,np) = [t trans].';
    else
        transitions(:,n_transitions_per_file) = [t trans].';
        %res.(['transitions' num2str(n/n_transitions_per_file)]) = transitions;
        eval(['trans' num2str(n/n_transitions_per_file) ' = transitions;'])
        save(transitionsFileName, ['trans' num2str(n/n_transitions_per_file)], '-append')
        disp(['Saving batch' num2str(n/n_transitions_per_file)])
        %inp = load(transitionsFileName, 'res');
        %pause()
    end
    n = n + 1;
end
%catch
%    finalState = state;
%    save(transitionsFileName, 'finalState', '-append')
    
%end
finalState = state;   
save(transitionsFileName, 'finalState', '-append')
