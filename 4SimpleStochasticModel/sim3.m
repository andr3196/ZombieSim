function sim3()
%% Goal 
% Make histogram animation for the enclave size of human survivers 
close all

fromFile = true;
file = 'data/enclavesScan2Power2Edges.mat';
nAlphaValues = 30 + 1;
alphaValues = linspace(0, 0.54, nAlphaValues);
edges = [0 2.^(0:8)];

if fromFile
    inp = load(file);
    enclaves = inp.enclaves;
else
    enclaves = getEnclaves(nAlphaValues,alphaValues, edges);
    save(file, 'enclaves')
end

% ignore alpha < 0.1
filt = alphaValues > 0.1;



bar3(alphaValues(filt),enclaves(filt, 2:end-1))
xlabel('N enclaves')
ylabel('\alpha')
zlabel('frequency')
view([138 19])
axis([1 9 0.1 0.54 0 0.8])
xticklabels({'1', '2', '4', '8', '16', '32', '64', '128'})
xtickangle(-45)
grid on
return
for i = 1: nAlphaValues
    bar(edges(1:end-1) + 5, enclaves(i,:))
    title(['\alpha = ' num2str(alphaValues(i))])
    axis([0 200 0 60])
    pause()
end


function enclaves = getEnclaves(nAlphaValues, alphaValues, edges)
    % grid size
n = 100;

%% Define scan
nReps = 50;
enclaves = zeros(nAlphaValues,length(edges)-1);

encCollect = zeros(nReps,length(edges)-1);

for j = 1:nAlphaValues
    alpha = alphaValues(j);
    
    
    for k = 1:nReps
        M = gridPropagate(n, [50,50], alpha);
        %% Calculate fractions
        SusMatrix = M == 0;
        CCSus = bwconncomp(SusMatrix);
        EnclaveSize = zeros(1,CCSus.NumObjects);
        for i = 1:CCSus.NumObjects
            EnclaveSize(i) = numel(CCSus.PixelIdxList{i});
        end 
        encCollect(k,:) = histcounts(EnclaveSize, edges, 'Normalization', 'probability');
    end
    enclaves(j,:) = mean(encCollect,1);
    
end










