function updateFigure(state, handles, has_changed)
%% Goal
% Updates the handles of the zones corresponding to the columns of the
% state that has changed
colorMat = getColorMat(state(1:3,:));

for i = 1:length(handles)
    if has_changed(i)
        set(handles(i), 'FaceColor', colorMat(:,i))
    end
end



function colorMat = getColorMat(state)
% Produces a color vector for each column in state
        colorMat = state./sum(state,1);
        colorMat = [0 1 0; 0 0 1;1 0 0]*colorMat; % blue = susc, green = zom, red = expo
        filter = ~isfinite(colorMat(1,:));
        n_filt = sum(filter);
        colorMat(:,filter) = repmat([1;1;1],1,n_filt); % make empty zones white 