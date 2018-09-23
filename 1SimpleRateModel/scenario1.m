function scenario1()
%% Goal 
% Perform a integration of the simple rate differential equation.

figure

N = 270000; % total population
Z0 = 1; % Initial zombie population

state = [N-Z0,Z0, 0].';

[t,state] = propagate(state);


plot(t,state(:,1), 'b-',t,state(:,2), 'g-')
legend({'Susceptibles', 'Zombies'})
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',extractValues('propagate.m'),'FitBoxToText','on')
xlabel('Time/hr')
ylabel('Number of humans/zombies')




