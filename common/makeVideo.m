function makeVideo(frames, filename, framerate)

myVideo = VideoWriter(filename, 'MPEG-4');
myVideo.FrameRate = framerate;  % Default 30
open(myVideo)
writeVideo(myVideo, frames);
close(myVideo)

