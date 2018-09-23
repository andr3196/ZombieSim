

close all
inp = load('run1scenario3M.mat');
M = inp.M;
k = 5;
n = 20;
theta = pi*(-n:2:n)/n;
phi = (pi/2)*(-n:2:n)'/n;

U = cos(phi)*cos(theta);
V = cos(phi)*sin(theta);
W = sin(phi)*ones(size(theta));
X = zeros(size(U));
Y = X;
Z = X;

hold on
for i = 1:size(U,1)
    for j = 1:size(U,2)
        if U(i,j)>=0 && V(i,j) >= 0 && W(i,j) >= 0
            C = [0 1 0; 0 0 1;1 0 0]* [U(i,j); V(i,j); W(i,j)];
            quiver3(X(i,j), Y(i,j), Z(i,j), U(i,j), V(i,j), W(i,j), 'Color', C)
        end
    end
end
grid on
view([168 22])
xlabel('Susceptibles')
ylabel('Infected')
zlabel('Zombies')