function batches = fetch_batches(transitions_file, batches_file, timeOffset)
%% Goal 
% Determine whether to recalculate batches or load from file

if ~isdefined('timeOffset', who)
    timeOffset = 0;
end


if ~exist(batches_file, 'file')
    batches = get_batches(transitions_file,timeOffset);
    save(batches_file, 'batches')
else
    inp = load(batches_file);
    batches = inp.batches;
end