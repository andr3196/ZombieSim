function zombiesInDenmarkPoly()
close all

text = fileread('DenmarkMap.rtf');

curStart = 1;
curEnd = NaN;

iMax = 5000;
i = 0;

figure
hold on
axis equal
n = 1;

% data: 
%name
%color
%polygon

data = cell(3,2160);

while i < iMax
    
    % find new path
    curStart = regexp(text, '<path ', 'once');
    if isempty(curStart)
        break
    end
    curEnd = regexp(text, '</path>', 'once');
    pathStr = text(curStart:curEnd);
    
    try
    % find color
    colorVecStart = regexp(pathStr, 'rgb(', 'once');
    colorVecEnd = regexp(pathStr, ')', 'once');
    colorVec = str2num(pathStr(colorVecStart + 4 :colorVecEnd - 1))/255;
    
    % find name
    nameTagStart = regexp(pathStr,'highcharts-name-');
    nameStart = nameTagStart + length('highcharts-name-');
    name = pathStr(nameStart:end-3);
    name = strrep(strrep(strrep(name, '\''f8','ø'), '\''e6', 'æ'), '\''e5', 'å');
    
    % find path values
    i1 = regexp(pathStr,'M', 'once');
    i2 = regexp(pathStr,'Z', 'once');
    pathValuesStr = pathStr(i1+1:i2-1);
    pathValuesStr(pathValuesStr == 'L') = '';
    pathValuesStr(pathValuesStr == 'M') = '';
    pathValues = str2num(pathValuesStr);
    %pathValues = [pathValues pathValues(1:2)];
    pathValues = [pathValues(1:2:end).' pathValues(2:2:end).'];
    pathValues = unique(pathValues, 'rows', 'stable');
    %pathValues = [pathValues; pathValues(1:2)];
    p = polyshape(pathValues);
    data(:,n) = {name, colorVec, p}.'; 
    n = n + 1;
    catch 
        %disp(['Could not parse: ' pathValuesStr])
        %if isempty(pathValuesStr)
         %   disp(pathStr)
        %end
    end 
    text = text(curEnd + 6:end);
    i = i + 1;
end
%set(gca,'Ydir','reverse')
save('DenmarkMap.mat', 'data')



% <path fill="rgb(99,169,218)" d="M 758.0063761811549 213.31877560775544 L 757.8144932428464 213.77612967770088 757.5120518148025 213.90717155300442 757.7970915104944 214.3521145475804 757.3673973046969 214.79626034460188 756.786252002531 214.09854214213772 757.3948166302122 213.5200684626846 757.742189937732 213.20753830024037 758.0063761811549 213.31877560775544 Z" stroke="#FFFFFF" stroke-width="0.5" vector-effect="non-scaling-stroke" class="highcharts-point highcharts-color-0 highcharts-name-helligånds">