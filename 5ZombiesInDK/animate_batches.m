function [state, frames] = animate_batches(state, batches, zoneHandles, pause_func, should_film)

%% Batch animate

frames(length(batches)) = struct('cdata',[],'colormap',[]);

for j = 1:length(batches)
    batch = batches{j};
    if ~isempty(batch)
        title(['Current time: ' num2str(batch{1}) ' hr'])
        [state, has_changed] = apply_batch(state, batch);
        updateFigure(state, zoneHandles, has_changed);
        pause_func()
        if should_film; frames(j) = getframe(gcf); end
    end
end