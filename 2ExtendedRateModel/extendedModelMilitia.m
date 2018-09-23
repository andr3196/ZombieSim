function extendedModelMilitia()
figure
global beta sigma betap sigmap kappa kappap rho delta;

% Base interaction rate

%%% PARAM BEGIN %%%
r = 0.01/91; % per hr
q = 0.8;
f = 10; %
beta = r/(2+q);
sigma = r/(2+q);
kappa = q*beta; 
lambda = 0;
nu = 0;
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

runningTime = [0 100];
nTSteps = 500;

% Setup trigger
opts = odeset('Events', @eventTrigger);

% perform scan
nu_values = 0;%.3:5e-3:0.6;
lambda_values = 0;%1e-7:1e-6:1e-4;

% Collect
scans = zeros(nTSteps, 5,length(nu_values));

for i = 1:length(lambda_values)
        nu = nu_values;
        lambda = lambda_values(i);
        state = [N-Z0-M0,0,Z0, M0, 0].';

        [t,state] = ode45(@(t,y) derivative(t,y, lambda, nu, tau),runningTime,state, opts);
        [t,state] = equidistantializer(t,state, nTSteps);
        scans(:,:,i)  = state;
end 
pause()

p = plot(t, scans(:,1,1) + scans(:,4,1) + scans(:,5,1), 'b-',t, scans(:,2,1), 'y-', t, scans(:,3,1),'g-');
title('Forming a militia')
l = legend({'Alive', 'Infected', 'Zombies'});
xlabel('Time/hr')
ylabel('Number of humans/zombies')
rh=annotation('textbox',l.Position - [0.0 0.2 0 0], 'String',endStatusString(sum(scans(end,[1 4 5],1)),scans(end,3,1),scans(end,4,1)),'FitBoxToText','on'); 

%frames(1) = getframe(gcf); 
for i = 2: length(lambda_values)
    title(['\lambda = ' num2str(lambda_values(i))])
    rh.String = endStatusString(sum(scans(end,[1 4 5],i)),scans(end,3,i));
    set(p,{'YData'}, num2cell(scans(:,:,i).',2))
    pause(0.25)
    %plot(t, scans(:,1,i), 'b-',t, scans(:,2,i), 'y-', t, scans(:,3,i),'g-', t, scans(:,4,i),'k-', t, scans(:,5,i),'r-')
    %frames(i) = getframe(gcf);
end
return 
myVideo = VideoWriter('extendedModelScanNuVid', 'MPEG-4');
myVideo.FrameRate = length(nu_values)/10;  % whole vid should take 10 s
%my%Video.Quality = 100;
open(myVideo)
writeVideo(myVideo, frames);
close(myVideo)
        


function dpdt = derivative(t, state, lambda, nu, tau)
global beta sigma sigmap betap rho delta kappap kappa
S = state(1);
I = state(2);
Z = state(3);
M = state(4);
V = state(5);

dSdt = (-(beta + sigma)*Z - nu + lambda*I)*S - tau*M*(1-M/(0.66*S+M));
dIdt = (beta*S + betap*M - sigma*I)*Z - lambda*S*I - rho*I;
dZdt = (-kappa*S - delta - kappap*M)*Z + rho*I;
dMdt = ( - sigmap*Z - nu - betap*Z)*M + tau*M*(1-M/(0.66*S+M));
dVdt = nu*(M+S) - sigma*V*Z;

dpdt = [dSdt; dIdt; dZdt; dMdt; dVdt];

 
function [value,isterminal,direction] = eventTrigger(~,state)
    value = ((state(1) + state(4) + state(5))-1);
    isterminal = 1;
    direction = 0;
    

    function s = endStatusString(nS, nZ, nMil)
        s = {'In the end:', ['Alive: ' num2str(floor(nS))], ['Zombies: ' num2str(floor(nZ))], ['Military: ' num2str(floor(nMil))] };
    
%% Runs


% Too fast
% PI = 0; % birth rate
% beta = 0.0095; % rate of Human -> Zombie
% ki = 0.0001; % Ressurection rate
% delta = 0.0001; % natural death rate
% alpha = 0.005; % rate Zombie -> Removed


