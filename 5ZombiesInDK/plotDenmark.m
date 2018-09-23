function [fig, handles] = plotDenmark()
%% Goal
%Plots map of Denmark
%% Load data
inpMap = load('data/DenmarkMapWithInhab.mat');
zoneData = inpMap.data;

%% Plot

f = figure;
    f.Name = 'ZOMBIESINDENMARK';
    hold on
    handles = gobjects(1,size(zoneData,2));
    for k = 1:size(zoneData,2)
        handles(k) = plot(zoneData{3,k}, 'FaceColor', 'none');
        
    end
    set(gca,'Ydir','reverse')

    fig = gcf;