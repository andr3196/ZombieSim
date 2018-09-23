function zombiesInDenmark5()
inp = load('zomInDKrun1.mat', 'res');
transitionStruct = inp.res;


%% Initialize state
inpMap = load('DenmarkMapWithInhab.mat');
zoneData = inpMap.data;

numInhab0 = [zoneData{4,:}];
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

%% Plot Denmark
curFig = gcf;
if ~strcmp(curFig.Name, 'ZOMBIESINDENMARK')
    f = figure;
    f.Name = 'ZOMBIESINDENMARK';
    hold on
    zonePlotHandles = gobjects(1,size(zoneData,2));
    for k = 1:size(zoneData,2)
        zonePlotHandles(k) = plot(zoneData{3,k}, 'FaceColor', 'none');
    end
    set(gca,'Ydir','reverse')
    savefig('CurrentZomInDKPlot')
    save('zombiePlotHandlesRun1.mat','zonePlotHandles');
else
    inpHandles = load('zombiePlotHandlesRun1.mat','zonePlotHandles');
    zonePlotHandles = inpHandles.zonePlotHandles;
    openfig('CurrentZomInDKPlot', 'visible')
end

%% Initialize zombie graphics



%% Prepare animation
fields = fieldnames(transitionStruct);
transitionNames = fields(contains(fields, 'transition'));


%% Animate
for i = 1:length(transitionNames)
    transitionName = transitionNames{i};
    transitions = transitionStruct.(transitionName);
    disp(transitionName)
    for j = 1:size(transitions,2)
        trans = transitions(:,j);
        title(['Current time: ' num2str(trans(1)) ' hr'])
        switch trans(2)
            case 1 % bite
                state(1:2,trans(3)) = state(1:2,trans(3)) + [-1; 1];
                colorVec = getColorVec(state(:,trans(3)));
                set(zonePlotHandles(trans(3)), 'FaceColor', colorVec)
            case 2 % transform
                state(2:3,trans(3)) = state(2:3,trans(3)) + [-1; 1];
                colorVec = getColorVec(state(:,trans(3)));
                set(zonePlotHandles(trans(3)), 'FaceColor', colorVec)
            case 3 % kill
                state(3:4,trans(3)) = state(3:4,trans(3)) + [-1; 1];
                colorVec = getColorVec(state(:,trans(3)));
                set(zonePlotHandles(trans(3)), 'FaceColor', colorVec)  
            case 4 % move
                state(3,[trans(3) trans(4)]) = state(3,[trans(3) trans(4)]) + [-1 1];
                colorVec1 = getColorVec(state(:,trans(3)));
                colorVec2 = getColorVec(state(:,trans(4)));
                set(zonePlotHandles(trans(3)), 'FaceColor', colorVec1) 
                set(zonePlotHandles(trans(4)), 'FaceColor', colorVec2)
        end
        
        pause(0.01)
    end
end


function colorVec = getColorVec(zoneState)

colorVec = zoneState(1:3)/sum(zoneState(1:3));
colorVec = [0 1 0; 0 0 1;1 0 0]*colorVec;
if any(~isfinite(colorVec))
   colorVec = [1 1 0];
end






