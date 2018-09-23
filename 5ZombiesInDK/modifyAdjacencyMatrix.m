function modifyAdjacencyMatrix()
%% Add bridge connections to adjacency matrix

%% Load old
inp = load('data/adjacencyMatrix.mat');
spar_adjacency_matrix = inp.spar_adjacency_matrix;

%% Define connections (zoneIndexFrom, zoneIndexTo)

% Lillebælt
%736 -> 2020

% Storebælt
% 638 -> 1929

% Limfjorden
% 1251 -> 1256
% 1250 -> 1255

bridge_connections = [736, 2020;
    638,1929;
    1251, 1256;
    1250, 1255];

for i = 1: size(bridge_connections,1)
    spar_adjacency_matrix(bridge_connections(i,1), bridge_connections(i,2)) = 1;
end

save('data/adjacencyMatrix.mat', 'spar_adjacency_matrix')
    


