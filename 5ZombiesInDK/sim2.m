function sim2()
%% Goal
% Here we batch animate a sim of zombies in DK

%% Setup
should_plot = true;
pause_length = 'none';

if isnumeric(pause_length)
    pause_func = @() pause(pause_length);
else
    pause_func = @() drawnow;
end

figure
%% Fetch data
transitions_file1 = 'data/zomInDKrun3.mat';
batches_file1 = 'data/zomInDKrun3Hour1Batches.mat';

batches = fetch_batches(transitions_file1, batches_file1);

%% Initialize
inp = load('data/zomInDKrun1.mat', 'res');
res = inp.res;
state = res.initState;

%% Initialize plot
if should_plot
    [~, zoneHandles] = plotDenmark();
    updateFigure(state, zoneHandles, ones(length(zoneHandles)))
end
pause()

%% Animate
[state, frames1] = animate_batches(state, batches, zoneHandles, pause_func, true);
%[~, frames2] = animate_batches(state, batches2, zoneHandles, pause_func, true);

%% Make animation
frames = frames1;
makeVideo(frames, 'videos/zombiesInDKFullMovie', length(frames)/60)

%% Create score
score = sum(state,2);
sum(score)
annotation('textbox',[0.7 0.7 0.2 0.2], 'String', {'Status:', ['Susceptible: ' num2str(score(1)) ], ['Infected: ' num2str(score(2))], ['Zombies: ' num2str(score(3))], ['Removed: ' num2str(score(4))] }, 'FitBoxToText','on' )



