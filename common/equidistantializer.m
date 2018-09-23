function [tNew, stateNew] = equidistantializer(t,state, nSteps)
tNew = linspace(min(t),max(t),nSteps);
stateNew = spline(t,state.',tNew);
stateNew = stateNew.';