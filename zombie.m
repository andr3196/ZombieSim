function zombie(input)

%% The parameters of a zombie

% Do you have to get bitten to become a zombie?

% Appearance: Whether they can close to humans 
% Infection time: 
% Speed: walking/running
% Killed by: blow to the head/ any way/ must be burned
% Social: Groups/ individuals
% Intelligence: random/controlled behaviour
% Transferability:
%   -> transfer rate:
%   -> 


%% Prepare for new sim
close all
 
 
%% Define 
global dx dy xmin xmax ymin ymax xframemin xframemax yframemin yframemax
global X Y
global runningTime nTSteps pauseStep
dx = input.dx; %50; % sheet resolution in x-direction
dy = input.dy; %50; % sheet resolution in y-direction
% XY bounds of the crowd initially
xmin = input.xmin; %-5;
xmax = input.xmax; %5;
ymin = input.ymin; %-5;
ymax = input.ymax; %5;

% XY bounds of the frame
xframemin = input.xframemin;
xframemax = input.xframemax;
yframemin = input.yframemin;
yframemax = input.yframemax;

nx = linspace(xmin,xmax,dx);
ny = linspace(ymin,ymax,dy);
[X, Y] = meshgrid(nx,ny);


%% Parameters

numZom = 1000;
zombieCenter = [500,500];

s = 50;
Z = numZom/(2*pi*s^2)*exp(-(((X - zombieCenter(1)).^2+(Y-zombieCenter(2)).^2)/(2*s^2)));



 
%% Simulation constants
pauseStep = input.pauseStep;% 0.05;
runningTime  = input.runningTime;% [0 50];
nTSteps = input.nTSteps; % 500;
 
%% PARAMETERS
 
 
% Let the state vector be defined such that
state(:,:,1) = Z;
%state(:,:,2) = dZdt
 
% Choose sim
run = 1;
switch run
    case 1
        simulate(state,input);  
    case 2
        scenario1()
end
 
 
 
 
 
 
function simulate(state, input) 


% Define time of calculation
global nTSteps runningTime
[t,state] = ode45(@(t,y) derivative(t,y),runningTime,toColumn(state));
[t,state] = equidistantializer(t,state, nTSteps);
animate(t,state)

 
function dpdt = derivative(~, state)
state = toColumnInverseN(state,1);
D = 10;
Z = state(:,:,1);
dZdt = D*del2(Z); 
dpdt = cat(3,dZdt);
dpdt = toColumn(dpdt);
 
 
%%% Helper functions
 
% Assumes M has dimensions n x m x 4
%function MAsColumn = toColumn(M)
%MAsColumn = reshape(M,numel(M),1,1);
 
%function  MAsMatrix = toColumnInverse(M)
%global dx dy
%MAsMatrix = reshape(M,dy, dx, 4);
 

%function  MAsMatrix = toColumnInverseN(M,n)
%global dx dy
%MAsMatrix = reshape(M,dy, dx, n);