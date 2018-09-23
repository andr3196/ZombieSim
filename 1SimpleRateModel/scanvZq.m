function scanvZq()
%% Goal 
% Perform scan of zombie velocity and q-factor. Plot as a 2-D surf
vZ = 0.1:0.1:12;
q = 0.9:0.001:1;

%% Grid
[XvZ, Yq] = meshgrid(vZ, q);
Z = zeros(size(XvZ));

%% Intital state
N = 270000; % total population
Z0 = 1; % Initial zombie population
initState = [N-Z0,Z0, 0].';

for i = 1:length(q)
    for j = 1:length(vZ)
        [~, state] = propagate(initState,'q',q(i),'vZ', vZ(j), 'runningTime', [0 24]);
        Z(i,j) = state(end,1);
    end
end

surf(XvZ, Yq, Z/270000, 'EdgeColor', 'none')
view([0 90])
xlabel('Zombie speed km/hr')
ylabel('Strength factor, q')
cb = colorbar;
cb.Label.String = 'S_{1 day}/S_0';
    

