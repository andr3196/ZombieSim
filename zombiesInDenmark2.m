close all

inp = load('DenmarkMap.mat', 'data');
[parishNInhab, txt] = xlsread('parishInhabitants2017.xlsx');
data = inp.data;
data{4,2180} = {};

parishNames = txt(4:end,2);

%% Find color to inhabitants mapping

for i = 1:2160
    if ~isempty(data{4,i})
        continue
    end
    dataName = data{1,i};
    found = false;
    for j = 1:length(parishNames)
        index_of_paren = strfind(parishNames{j}, '(');
        parishName = parishNames{j}(1:index_of_paren-1);
        dataNameVariations = makeAllVariations(dataName);
        
        for k = 1:length(dataNameVariations)
            if contains(parishName, dataNameVariations{k}, 'IgnoreCase', true)
                found = true;
                break
            end
        end
        if found
            break
        end
    end
    
    if found
        data{4,i} = parishNInhab(j);
        parishNInhab(j) = [];
        parishNames(j) = [];
    else
        data{4,i} = data{4,i-1};
    end
       
end



%save('DenmarkMapWithInhab.mat', 'data')
    
% t -th
% % - s
% u - v
% aa - å
% - - ' '

function variations = makeAllVariations(name)
variations = {name};
if contains(name, '-')
    variations = [variations strrep(variations, '-', ' ')];
end
if name(end) == 's'
   var2 = variations;
   for i = 1:length(var2)
       var2{i} = var2{i}(1:end-1);
   end
   variations = [variations var2];
end
if contains(name, 'aa')
    variations = [variations strrep(variations, 'aa', 'å')];
elseif contains(name, 'å')
    variations = [variations strrep(variations, 'å', 'aa')];
end
if contains(name, 'v')
    variations = [variations strrep(variations, 'v', 'u')];
end
end


% colors = [0.776470588235294 0.866666666666667 0.956862745098039;
%     0.768627450980392 0.862745098039216 0.952941176470588;
%     0.807843137254902 0.886274509803922 0.964705882352941;
%     0.580392156862745 0.764705882352941 0.905882352941177;
%     0.301960784313725 0.615686274509804 0.831372549019608;
%     0.435294117647059 0.686274509803922 0.866666666666667;
%     0.0705882352941177 0.494117647058824 0.772549019607843;
%     0.760784313725490 0.858823529411765 0.952941176470588;
%     0.454901960784314 0.694117647058824 0.870588235294118;
%     0 0.454901960784314 0.752941176470588
%     ];
% 
% inhabitants = [441; 522; 45; 2911; 6289;4678;9025; 12905; 4434; 15927];
% 
% plot( colors(:,1), inhabitants, 'k.')



%figure 
%hold on
%for i = 1:2159
%    c = data{2,i};
%    plot3(c(1),c(2), c(3),'k.');% c, 'Marker', '.')
    
%end