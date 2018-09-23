function animate(times,state)
global xframemin xframemax yframemin yframemax pauseStep X Y
s = toColumnInverseN(state(1,:),1);
h = surf(X,Y,s(:,:,1), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
axis([xframemin xframemax yframemin yframemax -0.01 1])
view([12 40])
%caxis([0 100]) 
colorbar
xlabel('x')
ylabel('y')
pause()
for i = 2:length(times)
    s = toColumnInverseN(state(i,:),1);
    set(h,'ZData', s(:,:,1));
    pause(pauseStep)
end