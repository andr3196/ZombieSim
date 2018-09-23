function extendedModelMilitiaVaccine()
figure
global beta sigma betap sigmap kappa kappap rho delta t0;

% Base interaction rate

%%% PARAM BEGIN %%%
t0 = 0; % Time at which we start spreading the vaccine
r = 0.01/91; % per hr
q = 0.8;
f = 10;
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

runningTime = [0 130];
nTSteps = 1000;

% Setup trigger
opts = odeset('Events', @eventTrigger);

% perform scan
nu_values = 0:1e-3:0.015;
lambda_values = (0:1e-3:0.015)*10; %0:1e-6:1e-4;

% Collect
scansNu = zeros(nTSteps, 5,length(nu_values));
scansLambda = zeros(nTSteps, 5,length(lambda_values));

fig.folder = 'plots';
fig.figuresModified = 12;
fig.filetype = '-dpdf';
fig.slimOff = true;


% Perform scan over vaccination rates
for i = 1:length(nu_values)*0
        nu = nu_values(i);
        lambda = 0;
        state = [N-Z0-M0,0,Z0, M0, 0].';

        [t,state] = ode45(@(t,y) derivative(t,y, lambda, nu, tau),runningTime,state, opts);
        [t,state] = equidistantializer(t,state, nTSteps);
        scansNu(:,:,i)  = state;
        i
        figure
        %set(gcf, 'Position', [500 500 600 500])
        p = plot(t, scansNu(:,1,i) + scansNu(:,4,i) + scansNu(:,5,i), 'b-',t, scansNu(:,2,i), 'y-', t, scansNu(:,3,i),'g-',t, scansNu(:,5,i),'r-');
        axis([runningTime 0 275e3])
        title(['\nu = ' num2str(nu*100) '%/hr']);
        l = legend({'Alive', 'Infected', 'Zombies', 'Vaccinated'});
        xlabel('Time/hr')
        ylabel('Number of humans/zombies')
        rh=annotation('textbox',l.Position - [0.03 0.2 0 0], 'String',endStatusString(sum(scansNu(end,[1 4 5],i)),scansNu(end,3,i),scansNu(end,4,i), scansNu(end,5,i)),'FitBoxToText','on'); 
        fig.filename = ['vaccine' num2str(i)];
        errorBarPlot(fig)
end


titles = {'No cure', 'With cure'};
% Perform scan over cure rates
for i = 1:2
        nu = 0;
        lambda = i - 1;
        state = [N-Z0-M0,0,Z0, M0, 0].';

        [t,state] = ode45(@(t,y) derivative(t,y, lambda, nu, tau),runningTime,state, opts);
        [t,state] = equidistantializer(t,state, nTSteps);
        scansLambda(:,:,i)  = state;
        i
        figure
        %set(gcf, 'Position', [500 500 600 500])
        p = plot(t, scansLambda(:,1,i) + scansLambda(:,4,i) + scansLambda(:,5,i), 'b-',t, scansLambda(:,2,i), 'y-', t, scansLambda(:,3,i),'g-');
        axis([runningTime 0 275e3])
        
        title(titles{i});
        l = legend({'Alive', 'Infected', 'Zombies'});
        xlabel('Time/hr')
        ylabel('Number of humans/zombies')
        rh=annotation('textbox',l.Position - [0.03 0.2 0 0], 'String',endStatusString(sum(scansLambda(end,[1 4 5],i)),scansLambda(end,3,i),scansLambda(end,4,i), scansLambda(end,5,i)),'FitBoxToText','on'); 
        fig.filename = ['cure' num2str(i)];
        errorBarPlot(fig)
end 

return

title('Starting treatment after 0 hrs')

% Plot effect of vaccination on the left 
s1 = subplot(1,2,1);

p = plot(t, scansNu(:,1,1) + scansNu(:,4,1) + scansNu(:,5,1), 'b-',t, scansNu(:,2,1), 'y-', t, scansNu(:,3,1),'g-',t, scansNu(:,5,1),'r-');

l = legend({'Alive', 'Infected', 'Zombies'});
xlabel('Time/hr')
ylabel('Number of humans/zombies')
rh=annotation('textbox',l.Position - [0.0 0.1 0 0], 'String',endStatusString(sum(scansNu(end,[1 4 5],1)),scansNu(end,3,1),scansNu(end,4,1)),'FitBoxToText','on'); 

% Plot effect of cure on the right
s2 = subplot(1,2,2);
p2 = plot(t, scansLambda(:,1,1) + scansLambda(:,4,1), 'b-',t, scansLambda(:,2,1), 'y-', t, scansLambda(:,3,1),'g-');
l2 = legend({'Alive', 'Infected', 'Zombies'});
rh2=annotation('textbox',l2.Position - [0.0 0.1 0 0], 'String',endStatusString(sum(scansLambda(end,[1 4],1)),scansLambda(end,3,1),scansLambda(end,4,1)),'FitBoxToText','on');
pause()
frames(1) = getframe(gcf);

for i = 2: length(nu_values)
    % Update annotation 
    rh.String = endStatusString(sum(scansNu(end,[1 4 5],i)),scansNu(end,3,i),scansNu(end,4,i) );
    rh2.String = endStatusString(sum(scansLambda(end,[1 4 5],i)),scansLambda(end,3,i),scansLambda(end,4,i) );
    % Update values
    val = [sum(scansNu(:,[1 4 5],i),2), scansNu(:,2,i), scansNu(:,3,i), scansNu(:,5,i)];
    val2 = [sum(scansLambda(:,[1 4],i),2), scansLambda(:,2,i), scansLambda(:,3,i)];
    set(p,{'YData'}, num2cell(val.',2))
    set(p2,{'YData'}, num2cell(val2.',2))
    % Update titles
    title(s1,['\nu = ' num2str(nu_values(i)*100) ' %/hr'])
    %title(s2,['\lambda = ' num2str(lambda_values(i)*100) ' %/hr'])
    % Fix axes
    axis(s1,[runningTime 0 275e3])
    axis(s2,[runningTime 0 275e3])
    drawnow
    frames(i) = getframe(gcf);
end
 
myVideo = VideoWriter('videos/vaccineCureScan', 'MPEG-4');
myVideo.FrameRate = 1;
%my%Video.Quality = 100;
open(myVideo)
writeVideo(myVideo, frames);
close(myVideo)
        


function dpdt = derivative(t, state, lambda, nu, tau)
global beta sigma sigmap betap rho delta kappap kappa t0 
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
        s = {'In the end:', ['Alive: ' num2str(floor(nS))], ['Zombies: ' num2str(floor(nZ))], ['Military: ' num2str(floor(nMil))], ['Vaccinated: ' num2str(floor(nVac))] };
    
%% Runs


% Too fast
% PI = 0; % birth rate
% beta = 0.0095; % rate of Human -> Zombie
% ki = 0.0001; % Ressurection rate
% delta = 0.0001; % natural death rate
% alpha = 0.005; % rate Zombie -> Removed


