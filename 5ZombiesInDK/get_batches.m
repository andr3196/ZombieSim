function batches = get_batches(filename, timeOffset)
%% Goal
% Given a filename to load transitions data from, we create batches for every hour, batch_size, a batch being defined below
%% Load data
inp = load(filename);
names = fieldnames(inp);

%% Setup loop over 'trans*' fields
trans_tag = 'trans';

n_trans_segments = sum(contains(names,trans_tag)); % subtracting 1 'res' field


% To get the total number of batches with 1 batch/hr, we need the time for
% the last event in the last trans

% Get last trans
last_trans = inp.(['trans' num2str(n_trans_segments)]);

% The total number of hours is then (also num_batches
num_batches = ceil(last_trans(1,end));

%%%%num_trans = size(inp.trans1,2);
%%%%%num_batches = floor(num_trans/batch_size);

batches = cell(1,num_batches);

% A 'batch' has the format {t_start, t_end, case1, case2,..., case4b},
% where
% casei = [k1 k2... ;n1 n2 ], where
% kj is a zone number, and nj is the change in one of [S_kj, E_kj, Z_kj,
% R_kj] given
% case1: nj is the number of susc    who have been bitten
% case2:                     exposed               transformed
% case3:                     zombies               killed
% case4a nj                  zombies               left zone kj
% case4b                     zombies               entered zone kj.
% case5                      military              trained
% case6                      zombies               killed by mili
% case7                      military              bitten       
% case8                      military              entered zone kj

nCases = 8;

% The index of the trans field
cur_trans_index = 2;
cur_trans = inp.(names{cur_trans_index});
cur_batch_index = 1;
cur_start_index = 1;

if ~isdefined('timeOffset',who)
    timeOffSet = 0;
end

for t_start = timeOffset:timeOffset + num_batches - 1
    %% Init new hour
    t_end = t_start + 1;
    batch_events = [];
    %% Resolve batch events in current hour
    while true
        sep_index = find(cur_trans(1,cur_start_index:end) > t_end, 1); % index of first event not happening in this hour
        if sep_index == 1 % all events in the previous trans happen in this hour, but none in this one
            break
        elseif isempty(sep_index) % All events in cur_trans happen in this hour
            batch_events = [batch_events cur_trans];
            cur_trans_index = cur_trans_index + 1;
            if cur_trans_index <= n_trans_segments + 1
                cur_trans = inp.(names{cur_trans_index});
            else
                break
            end
            cur_start_index = 1;
        else % only some events in cur_trans happen in this hour
            batch_events = [batch_events cur_trans(:, cur_start_index:cur_start_index + sep_index - 2 ) ];
            cur_start_index = cur_start_index + sep_index - 1;
            break
        end
    end
    %% Make batch
    % sort into batches
    [B,I] = sort(batch_events(2,:));
    batch = batch_events(:,I);
    % in the sequence "1 1 1 2 2 3 3 3 3 4..." find the indices where 1->2, 2->3... etc.
    type_separation_indices = find(diff(B));
    types = B(type_separation_indices);
    nTypes = length(type_separation_indices);
    
    batch_t = cell(1, 2 + nCases + 2);
    batch_t{1} = t_start;
    batch_t{2} = t_end;
    batch_t{3:2+nCases} = [];
    
    % events batch(ind1:ind2) belong 
    ind1 = 1;
    ind2 = type_separation_indices(1);
    for k = nTypes
        typ = types(k);
            if typ == 4 || typ == 8
                batch{2 + typ} = countOccurances(batch(3, ind1:ind2));
                batch{2 + typ + 1} = countOccurances(batch(4, ind1:ind2));
            else
            batch{2 + types(k)} = countOccurances(batch(3,ind1:ind2));
            end
            
            ind1 = ind2 + 1;
            
            ind2 = type_separation_indices(1)
            
    end
    
    
    
    if nTypes > 0
    %% case 1 - Bite
    batch_case1 = countOccurances(batch(3,1:type_separation_indices(1)));
    else 
        batch_case1 = [];
    end
    % transform from (i, + 1) to (i, +n)
    if nTypes > 1
    %% case 2 - Transform
    batch_case2 = countOccurances(batch(3, type_separation_indices(1)+1:type_separation_indices(2)));
    else
       batch_case2 = [];
    end
    
    if nTypes > 2
    %% case 3 - Kill
        batch_case3 = countOccurances(batch(3, type_separation_indices(2)+1:type_separation_indices(3)));
    else
        batch_case3 = [];
    end
    
    if nTypes > 3
    %% case 4 - Move
    batch_case4from = countOccurances(batch(3, type_separation_indices(3)+1:end));
    batch_case4to = countOccurances(batch(4, type_separation_indices(3)+1:end));
    else
        batch_case4from = [];
        batch_case4to = [];
    end 
    
    
    
    batches{cur_batch_index} = batch_t;
    
    cur_batch_index = cur_batch_index + 1;
    
end



function counts = countOccurances(list)
% in a list of zone numbers, kj, [1, 34, 1, 2, 3, 3,1...] count the number of
% occurances of each number and returns a 2xn matrix of [kj, ...; nj,...],
% where nj is the number of times kj occurs
zone_max = max(list);
counts = zeros(2,zone_max);
for zone = 1:zone_max
    counts(:,zone) = [zone, sum(list(:) == zone)];
end
counts(:, counts(2,:) == 0) = [];

