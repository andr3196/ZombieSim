function i = getParishByCoord(zoneData, p)

x = p(1);
y = p(2);

for i = 1:length(zoneData)
if isinterior(zoneData{3,i}, x, y)
break;
end
end