function [M, film, Ms] = gridPropagate(n, zomInitialPos, alpha, fixedRandom, shouldAnimate, shouldFilm, nthFrame, shouldCollectM)
%% Goal

%% Define
% Let
% 0 = Susceptible
% 1 = Zombie
% -1 = Removed
global n nBonds

%% Parameters
%beta = ; % zombie bite rate
%kappa = ; % human kill rate

vars = who;
if ~isdefined('alpha', vars)
    alpha = 0.4; %kappa/beta
end
if ~isdefined('n', vars)
    n = 100; % grid side length
end

if ~isdefined('zomInitialPos', vars)
    zomInitialPos = [1,1]; % Place of the first zombies
end

if ~isdefined('fixedRandom', vars)
    fixedRandom = false;
end

if ~isdefined('shouldAnimate', vars)
    shouldAnimate = false;
end

if ~isdefined('shouldFilm', vars)
    shouldFilm = false;
else
    if ~isdefined('nthFrame', vars)
    	nthFrame = 1; 
    end
end

if ~isdefined('shouldCollectM', vars)
    shouldCollectMs = false;
else
    sampleRate = 10;
    Ms = zeros(n,n,n^2/sampleRate);
    shouldCollectMs = true;
end



%% Setup


% Should seed random
if fixedRandom
    rng(22);
end

% Define grid
M = zeros(n);

% Initial zombie position (1,1)
M(zomInitialPos(1),zomInitialPos(2)) = 1;
M(n,n) = -1;

% Interaction queue: i1, j1 (zombie), i2, j2 (human)
nBonds = 0; % number of rows in Q
% Adjacent tiles
A = [0 1; 1 0; 0 -1; -1 0];
Q = appendBonds([],A,M,zomInitialPos); % nBonds x 4 matrix

% Count steps
step = 0;

if shouldAnimate
    % Plot
    map = [1 0 0;
        0 0 1;
        0 1 0];
    h = imagesc(M);
    colormap(map)
    colorbar('Ticks',[-1, 0, 1],...
        'TickLabels',{'Removed', 'Susceptible', 'Zombie'})
    xlabel('x')
    ylabel('y')
    axis equal
    set(gca,'FontSize',15)
    pause()
    film(ceil(n^2/nthFrame)) = struct('cdata',[],'colormap',[]);
    film(1) = getframe(gcf);
end
if shouldCollectMs 
    Ms(:,:,1) = M;
end

rands = rand(1,4*n^2);

while true
    % Select a random bond
    randIndex = randi(nBonds);
    bond = Q(randIndex,:);
    % Remove bond from Q
    Q(randIndex,:) = [];
    nBonds = nBonds - 1;
    if rands(step + 1) < 1/(1+alpha) % human is bitten
        M(bond(3), bond(4)) = 1;
        Q = appendBonds(Q,A,M,bond(3:4));
    else % zombie is killed
        M(bond(1), bond(2)) = -1;
        Q = removeBonds(Q,bond(1:2));
    end
    if shouldCollectMs && mod(step+2,sampleRate) == 0 
        Ms(:,:,(step+2)/sampleRate) = M;
    end
    if shouldAnimate
        set(h, 'CData', M)
        if shouldFilm && mod(step, nthFrame) == 0
            drawnow
            film(step/nthFrame + 2) = getframe(gcf);
        else
            %pause(0.01)
        end
    end
    if nBonds == 0
         if shouldFilm
             film(step + 3:end) = [];       
         end
        break
    end
    
    step = step + 1;
end


if ~isdefined('film', who)
    film = [];
end
if ~isdefined('Ms', who)
    Ms =[];
end
    


function point = periodicBoundary(point)
%% Wraps the point around the boundaries of the grid such that 0->n, n+1->1
global n
point(point == 0) = n;
point(point == n+1) = 1;

function Q = appendBonds(Q,A,M,point)
%%
global nBonds
for i = 1:4
    newPoint = periodicBoundary(point + A(i,:));
    if M(newPoint(1),newPoint(2)) == 0
        Q(end+1,:) = [point newPoint];
        nBonds = nBonds + 1;
    end
end

function Q = removeBonds(Q, point)
global nBonds
filter = all(Q(:,1:2) == point,2);
nBonds  = nBonds - sum(filter);
Q(filter,:) = [];



    
