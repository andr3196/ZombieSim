function phase_plot()
%% Goal
% Plot the direction of change in S and Z for different values of S, Z

%%% PARAM BEGIN %%%
setup.vZ = 1;
setup.q = 0.8;
r = 0.01*setup.vZ/91; % per hr
setup.beta = r/(2+setup.q);
setup.sigma = r/(2+setup.q);
setup.alpha = setup.q*setup.beta; 
%%% PARAM END %%%

%% Initialize populations
nPop = 270000;
nZom = 0:15000:nPop;
nHum = 0:15000:nPop;

[S,Z] = meshgrid(nHum, nZom);

%% Limits
% There cannot be more than 270000 individuals total

f = (S + Z) > nPop;

S(f) = NaN;
Z(f) = NaN;

[dS, dZ] = derivative(S,Z, setup);

quiver(S,Z, dS,dZ);
title('S-Z phase plot - vZ = 1 km/hr and q = 0.8')
xlabel('Number of susceptibles')
ylabel('Number of zombies')



function [dSdt, dZdt] = derivative(S,Z, setup)
dSdt = -(setup.beta + setup.sigma)*S.*Z;
dZdt = (setup.beta - setup.alpha)*S.*Z; 