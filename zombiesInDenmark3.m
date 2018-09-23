
inp = load('DenmarkMapWithInhab.mat', 'data');
data = inp.data;
nParish = size(data,2);
is_adjacent = zeros(nParish);
centroids = cell(1,nParish);
for i = 1:nParish
    [x, y] = centroid(data{3,i});
    centroids{i} = [x y];
end

adj_dist_sq = 20^2; 

for i = 1:nParish
    p1 = data{3,i};
    n1 = p1.NumRegions;
    c1 = centroids{i}; 
    for j = i+1:nParish
       if sum((c1 - centroids{j}).^2) > adj_dist_sq
           continue
       end
           
       p2 = data{3,j};
       u = union(p1,p2);
       is_adjacent(i,j) = u.NumRegions < n1 + p1.NumRegions;
    end
    
    if mod(i,100) == 0
       disp(i); 
    end
end

save('adjacencyMatrix.mat', 'is_adjacent')