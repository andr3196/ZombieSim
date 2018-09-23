function zombie()

input.dx =100;
input.dy = 100;
input.xmin = 0;
input.xmax = 100;
input.ymin = 0;
input.ymax = 100;

input.pauseStep = 0.05;
input.runningTime = [0 50];
input.nTSteps = 500;

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
global X Y deltax deltay
global runningTime nTSteps pauseStep
dx = input.dx; %50; % sheet resolution in x-direction
dy = input.dy; %50; % sheet resolution in y-direction
% XY bounds of the crowd initially
xmin = input.xmin; %-5;
xmax = input.xmax; %5;
ymin = input.ymin; %-5;
ymax = input.ymax; %5;

% XY bounds of the frame
xframemin = input.xmin;
xframemax = input.xmax;
yframemin = input.ymin;
yframemax = input.ymax;

nx = linspace(xmin,xmax,dx);
ny = linspace(ymin,ymax,dy);
[X, Y] = meshgrid(nx,ny);
deltax = nx(2)-nx(1);
deltay = ny(2)-ny(1);


%% Parameters

numZom = 10;
zombieCenter = [80,80];

sZ = 3;
Z = numZom/(2*pi*sZ^2)*exp(-(((X - zombieCenter(1)).^2+(Y-zombieCenter(2)).^2)/(2*sZ^2)));

numHuman = 100;
humanCenter = [70,70];
sH = 10;
H = numHuman/(2*pi*sH^2)*exp(-(((X - humanCenter(1)).^2+(Y-humanCenter(2)).^2)/(2*sH^2)));

 
%% Simulation constants
pauseStep = input.pauseStep;% 0.05;
runningTime  = input.runningTime;% [0 50];
nTSteps = input.nTSteps; % 500;
 
%% PARAMETERS
 
 
% Let the state vector be defined such that
state(:,:,1) = Z;
state(:,:,2) = H;
 
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
global deltax deltay
state = toColumnInverseN(state,2);
D = 20;
DH = 5;
beta = 3;
k = 0.5;
Z = state(:,:,1);
H = state(:,:,2);
dZdt = D*del2(Z, deltax, deltay) + beta*H.*Z - DH*del2(H, deltax, deltay); 
dHdt = -beta*H.*Z;
dpdt = cat(3,dZdt, dHdt);
dpdt = toColumn(dpdt);

function animate(times,state)
global xframemin xframemax yframemin yframemax pauseStep X Y
s = toColumnInverseN(state(1,:),2);
%ax1 = subplot(2,1,1);
%z = surf(X,Y,s(:,:,1), 1e3*s(:,:,1), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%colormap summer
h = imagesc(diff(s,1,3));
hold on
%ax2 = subplot(2,1,2);
%h = surf(X,Y,s(:,:,2), -1e3*s(:,:,2), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%colormap winter
%axis([xframemin xframemax yframemin yframemax -0.01 1])
%view([12 40])
colorbar
xlabel('x')
ylabel('y')
pause()
for i = 2:length(times)
    s = toColumnInverseN(state(i,:),2);
    %set(z,'ZData', s(:,:,1));
    set(h,'CData', diff(s,1,3));
    %pause(pauseStep)
    s 
    pause()
end


sum(sum(s(:,:,1)))*100^2
sum(sum(s(:,:,2)))*100^2

 
 
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