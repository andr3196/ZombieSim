function extendedModelCureTimeScan()
close all
global beta sigma betap sigmap kappa kappap rho delta;

% Base interaction rate

%%% PARAM BEGIN %%%
%t0 = 0; % Time at which we start spreading the vaccine
r = 0.01/91; % per hr
q = 0.8;
f = 10;
beta = r/(2+q);
sigma = r/(2+q);
kappa = q*beta; 
nu = 5e-3;
lambda = 0;
tau = 20/441; %/hr % One trained unit every month
betap = beta/f;
sigmap = sigma/f;
kappap = f*kappa;
delta = 1/(14*24); %/hr
rho = 2;% /hr
%%% PARAM END %%% 

N = 270000; % total population
Z0 = 1; % Initial zombie population
M0 = 10;

runningTime = [0 130];
nTSteps = 1000;

% Setup trigger
opts = odeset('Events', @eventTrigger);

% perform scan
times = [0 10.4 10.5 12];

% Collect
scansTimeVac = zeros(nTSteps, 5,length(times));
scansTimeCure = zeros(nTSteps, 5,length(times));

fig.folder = 'plots';
fig.filetype = '-dpdf';
fig.slimOff = 1;


% Perform scan over vaccination rates
for i = 1:length(times)*0
        state = [N-Z0-M0,0,Z0, M0, 0].';
        t0 = times(i);
        [t,state] = ode45(@(t,y) derivative(t,y,lambda, nu, tau, t0),runningTime,state);
        [t,state] = equidistantializer(t,state, nTSteps);
        scansTimeVac(:,:,i)  = state;
        i
        figure
        plot(t, scansTimeVac(:,1,i) + scansTimeVac(:,4,i) + scansTimeVac(:,5,i), 'b-',t, scansTimeVac(:,2,i), 'y-', t, scansTimeVac(:,3,i),'g-',t, scansTimeVac(:,5,i),'r-');
        axis([runningTime 0 275e3])
        title(['t_0 = ' num2str(times(i)) ' hr, \nu = 0.5 %/pr'])
        l = legend({'Alive', 'Infected', 'Zombies', 'Vaccinated'});
        xlabel('Time/hr')
        ylabel('Number of humans/zombies')
        rh=annotation('textbox',l.Position - [0.05 0.2 0 0], 'String',endStatusString(sum(scansTimeVac(end,[1 4 5],i)),scansTimeVac(end,3,i),scansTimeVac(end,4,i), scansTimeVac(end,5,i)),'FitBoxToText','on'); 
        fig.filename = ['vaccineTimeScan' num2str(i)];
        errorBarPlot(fig);
end
% Turn vaccine off, cure on
nu = 0;
lambda = 1;
% Perform scan over cure rates
for i = 1:length(times)
        t0 = times(i);
        state = [N-Z0-M0,0,Z0, M0, 0].';
        [t,state] = ode45(@(t,y) derivative(t,y,lambda, nu, tau, t0),runningTime,state);
        [t,state] = equidistantializer(t,state, nTSteps);
        scansTimeCure(:,:,i)  = state;
        i
        figure
        plot(t, scansTimeCure(:,1,i) + scansTimeCure(:,4,i) + scansTimeCure(:,5,i), 'b-',t, scansTimeCure(:,2,i), 'y-', t, scansTimeCure(:,3,i),'g-');
        axis([runningTime 0 275e3])
        title(['t_0 = ' num2str(times(i)) ' hr, with cure'])
        l = legend({'Alive', 'Infected', 'Zombies'});
        xlabel('Time/hr')
        ylabel('Number of humans/zombies')
        rh=annotation('textbox',l.Position - [0.05 0.2 0 0], 'String',endStatusString(sum(scansTimeCure(end,[1 4 5],i)),scansTimeCure(end,3,i),scansTimeCure(end,4,i), scansTimeCure(end,5,i)),'FitBoxToText','on'); 
        fig.filename = ['cureTimeScan' num2str(i)];
        errorBarPlot(fig);
        
end 
return
title('Starting treatment after 0 hrs')

% Plot effect of vaccination on the left 
s1 = subplot(1,2,1);

p = plot(t, scansTimeVac(:,1,1) + scansTimeVac(:,4,1) + scansTimeVac(:,5,1), 'b-',t, scansTimeVac(:,2,1), 'y-', t, scansTimeVac(:,3,1),'g-',t, scansTimeVac(:,5,1),'r-');
axis([runningTime 0 275e3])
title(s1,['t_0 = ' num2str(times(1)) ' hr'])
l = legend({'Alive', 'Infected', 'Zombies'});
xlabel('Time/hr')
ylabel('Number of humans/zombies')
rh=annotation('textbox',l.Position - [0.0 0.1 0 0], 'String',endStatusString(sum(scansTimeVac(end,[1 4 5],1)),scansTimeVac(end,3,1),scansTimeVac(end,4,1)),'FitBoxToText','on'); 

% Plot effect of cure on the right
s2 = subplot(1,2,2);
p2 = plot(t, scansTimeCure(:,1,1) + scansTimeCure(:,4,1), 'b-',t, scansTimeCure(:,2,1), 'y-', t, scansTimeCure(:,3,1),'g-');
l2 = legend({'Alive', 'Infected', 'Zombies'});
rh2=annotation('textbox',l2.Position - [0.0 0.1 0 0], 'String',endStatusString(sum(scansTimeCure(end,[1 4],1)),scansTimeCure(end,3,1),scansTimeCure(end,4,1)),'FitBoxToText','on');
axis([runningTime 0 275e3])
pause()
frames(1) = getframe(gcf);

for i = 2: length(times)
    % Update annotation 
    rh.String = endStatusString(sum(scansTimeVac(end,[1 4 5],i)),scansTimeVac(end,3,i),scansTimeVac(end,4,i) );
    rh2.String = endStatusString(sum(scansTimeCure(end,[1 4 5],i)),scansTimeCure(end,3,i),scansTimeCure(end,4,i) );
    % Update values
    val = [sum(scansTimeVac(:,[1 4 5],i),2), scansTimeVac(:,2,i), scansTimeVac(:,3,i), scansTimeVac(:,5,i)];
    val2 = [sum(scansTimeCure(:,[1 4],i),2), scansTimeCure(:,2,i), scansTimeCure(:,3,i)];
    set(p,{'YData'}, num2cell(val.',2))
    set(p2,{'YData'}, num2cell(val2.',2))
    % Update titles
    title(s1,['t_0 = ' num2str(times(i)) ' hr'])
    %title(s2,['\lambda = ' num2str(lambda_values(i)*100) ' %/hr'])
    % Fix axes
    axis(s1,[runningTime 0 275e3])
    axis(s2,[runningTime 0 275e3])
    drawnow
    frames(i) = getframe(gcf);
end
 
myVideo = VideoWriter('videos/cureTimeScan', 'MPEG-4');
myVideo.FrameRate = 1;
%my%Video.Quality = 100;
open(myVideo)
writeVideo(myVideo, frames);
close(myVideo)
close
        


function dpdt = derivative(t, state, lambda, nu, tau, t0)
global beta sigma sigmap betap rho delta kappap kappa
S = state(1);
I = state(2);
Z = state(3);
M = state(4);
V = state(5);
C = 1e-4*S;
dSdt = -(beta + sigma)*Z*S - (t>t0)*nu*S + (t>t0)*min([I,C])*(lambda>0) - tau*M*(1-M/(S+M));
dIdt = (beta*S + betap*M - sigma*I)*Z - (t>t0)*min([I,C])*(lambda>0) - rho*I;
dZdt = (-kappa*(S+V) - delta - kappap*M)*Z + rho*I;
dMdt = ( - sigmap*Z - betap*Z)*M + tau*M*(1-M/(S+M));
dVdt = (t>t0)*nu*S - sigma*V*Z;

dpdt = [dSdt; dIdt; dZdt; dMdt; dVdt];

 
function [value,isterminal,direction] = eventTrigger(~,state)
    value = ((state(1) + state(4) + state(5))-1);
    isterminal = 1;
    direction = 0;
    

function s = endStatusString(nS, nZ, nMil, nVac)
    nZ = max([nZ 0]);
    s = {'In the end:', ['Alive: ' num2str(floor(nS))], ['Zombies: ' num2str(floor(nZ))], ['Military: ' num2str(floor(nMil))], ['Vaccinated: ' num2str(floor(nVac))] };
    



