function [t, state] = propagate(initState, varargin)
%% Goal
% Propagate an initial state to the final time and produce a matrix equidistant points 

%% Setup
setup.runningTime = [0 10];
setup.nTSteps = 500;


state = initState;


opts = odeset('Events', @eventTrigger);

%%% PARAM BEGIN %%%
setup.vZ = 1;
setup.q = 0.8;
for i = 1:2:nargin-1
   setup.(varargin{i}) = varargin{i+1}; 
end
r = 0.01*setup.vZ/91; % per hr
setup.beta = r/(2+setup.q);
setup.sigma = r/(2+setup.q);
setup.alpha = setup.q*setup.beta; 
%%% PARAM END %%%

[t,state] = ode45(@(t,y) derivative(t,y, setup),setup.runningTime,state, opts);
[t,state] = equidistantializer(t,state, setup.nTSteps);


function dpdt = derivative(~, state, setup)
S = state(1);
Z = state(2);
dSdt = -(setup.beta + setup.sigma)*S*Z;
dZdt = (setup.beta - setup.alpha)*S*Z; 
dRdt = (setup.alpha + setup.sigma)*S*Z;
dpdt = [dSdt; dZdt; dRdt];

function [value,isterminal,direction] = eventTrigger(~,state)
    value = state(1)-1;
    isterminal = 1;
    direction = 0;

