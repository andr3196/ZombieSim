function scenario3()
%% Goal 
% Simulate stochastic spread of zombies on a grid


%% Prepare new figure 
figure

%% Run propagation
[M, frames] = gridPropagate([],[50, 50], [], true, true, true, 100);
%makeVideo(frames, 'videos/scen3run2', 30)














