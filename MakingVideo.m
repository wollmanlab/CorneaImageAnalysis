fpath = '/bigstore/Images/Alon/Analysis/EDOF0803/';
outputVideo = VideoWriter(sprintf('%smovie_1%s',fpath,'.avi'));
outputVideo.FrameRate = 1;
open(outputVideo)
%Loop through the image sequence, load each image, and then write it to the video.

%stk is a stack of images
for ii = 1:size(stk,3)
   writeVideo(outputVideo,stk(:,:,ii));
end
%Finalize the video file.

close(outputVideo)