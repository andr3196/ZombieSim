function inp = getInput()
%% Define input parameters
inp.dx = 50;
inp.dy = 50;
inp.xmin = 0;
inp.xmax = 1000;
inp.ymin = 0;
inp.ymax = 1000;
inp.xframemin = -10;
inp.xframemax = 1010;
inp.yframemin = -10;
inp.yframemax = 1010;
inp.pauseStep = 0.02;
inp.runningTime = [0 20];
inp.nTSteps = 500;