function findAvgParishSideLength()
close
%% Goal
% Assuming all Danish parishes are square - what is the average side
% length?
inp = load('data/DenmarkMapWithInhab.mat', 'data');
data = inp.data;
nParish = size(data,2);
areas = zeros(1,nParish);
for i = 1:nParish
    areas(i) = area(data{3,i});
end

% Real land area of Denmark
area_real = 42394; %km2

% normalize areas 
areas = areas*area_real/sum(areas);

% Find mean side length
sqrt(mean(areas))

% Plot histogram of areas
histogram(areas, 'Normalization', 'probability')
xlabel('Parish area/ km^2')
ylabel('Frequency')