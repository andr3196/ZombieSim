function sim2()
figure
%% Goal
% Plot the fraction of sites assigned zombie/susceptible as an average of
% 100 runs
close all

fromFile = true;
file = 'data/fractionsScanFullRange.mat';
nAlphaValues = 30 + 1;
alphaValues = linspace(0, 0.54, nAlphaValues);

if exist(file, 'file')
    inp = load(file);
    fractions = inp.fractions;
else
    fractions = getFractions(alphaValues, nAlphaValues);
    save(file, 'fractions')
end
%figure 
%for i = 1:nAlphaValues
%histogram(fractions(i,:,1), [0:0.1:1])
%pause
%end
%return

meanH = mean(fractions(:,:,1),2);
stdH = std(fractions(:,:,1),0,2);
meanZ = mean(fractions(:,:,2),2);
stdZ = std(fractions(:,:,2),0,2);

% pm 1 sigma human
Htop = meanH + stdH;
Hbottom = meanH - stdH;

% pm 1 sigma zombies 
Ztop = meanZ + stdZ;
Zbottom = meanZ - stdZ;


plot(alphaValues, meanH, 'bx')
hold on
fill ([alphaValues fliplr(alphaValues)], [Hbottom.' fliplr(Htop.')], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
%plot(alphaValues, Htop, 'r-')
%plot(alphaValues, Hbottom, 'r-')
plot(alphaValues,meanZ, 'Marker', '*','Color', [34, 139,34]/255, 'LineStyle', 'none')
fill ([alphaValues fliplr(alphaValues)], [Zbottom.' fliplr(Ztop.')], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
xlabel('Strength factor, q')
ylabel('Fraction of total population')
legend({'Susceptibles', '\pm \sigma (Susceptibles)' 'Zombies', '\pm \sigma (Zombies)'}, 'Location', 'north')
axis([0 0.54 0 1])





function fractions = getFractions(alphaValues,nAlphaValues)
    % grid size
n = 100;

%% Define scan
nReps = 100;

fractions = zeros(nAlphaValues, nReps, 2);


for j = 1:nAlphaValues
    alpha = alphaValues(j);
    
    for k = 1:nReps
        M = gridPropagate(n, [50,50], alpha);
        %% Calculate fractions
        fractions(j,k,1) = sum(sum(M == 0))/n^2;
        fractions(j,k,2) = sum(sum(M == 1))/n^2;  
    end
    
end








