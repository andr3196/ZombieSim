function [state, has_changed] = apply_batch(state, batch)
%% Goal
% changes the state according to the events listed in 'batch', returns also
% which zones have changed
oldState = state;
%% Case 1 - Bite: (Si, Ei) -> (Si - 1, Ei + 1)
bc1 = batch{3};
if ~isempty(bc1)
    state(1:2, bc1(1,:)) = state(1:2, bc1(1,:)) + [- bc1(2,:); bc1(2,:)];
end

%% Case 2 - Transform: (Ei, Zi) -> (Ei - 1, Zi + 1)
bc1 = batch{4};
if ~isempty(bc1)
    state(2:3, bc1(1,:)) = state(2:3, bc1(1,:)) + [- bc1(2,:); bc1(2,:)];
end
%% Case 3 - Kill: (Zi, Ri) -> (Zi -1, Ri +1)
bc1 = batch{5};
if ~isempty(bc1)
    state(3:4, bc1(1,:)) = state(3:4, bc1(1,:)) + [- bc1(2,:); bc1(2,:)];
end

%% Case 4 - Move: (Zi, Zj) -> (Zi -1, Zj +1)
bc1 = batch{6};
if ~isempty(bc1)
    state(3, bc1(1,:)) = state(3, bc1(1,:)) - bc1(2,:);
    bc1 = batch{7};
    state(3, bc1(1,:)) = state(3, bc1(1,:)) + bc1(2,:);
end

state(state < 0) = 0;

has_changed = any(oldState ~= state,1);