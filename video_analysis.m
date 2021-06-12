% Calculating the characteristics of a groundwater flow using a fixed camera
% Code finds the number on the ruler using template matching (normxcorr2 algorithm), references the frame accordingly
% Code finds the drop and evaluates its characterstics considering parallax correction 
% Written by: Mohammad Afzal Shadab
% Last edited: May 29th, 2021
% Email: mashadab@utexas.edu
% Input: raw video 
% Output: analysed video, datasheet
clc; % Clear the command window.
close all; % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear; % Erase all existing variables.
workspace; % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 14;
data = []; %the data will be compiled here
fps = 60; %frames rate (per second)
threshold = 70; %setting the threshold of the image


%Manually setting the region of importance
top_height = 321;
bottom_height = 543;
x_origin = 96;
x_right  = 1243;

crop_drop=  [x_origin top_height (x_right - x_origin) (bottom_height - top_height)]; %cropping the drop region: left, top, width, height
frame_begin = 790;  %starting frame
frame_end = 30790;   %end frame
frame_skip = 60;
position_timer = [0, 50];   %position of the timer

% Read the video in a standard MATLAB color video format
folder = fullfile('\');
fffilename = 'June8_1st_transient_try';
baseFileName = sprintf('%s.mov',fffilename);
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end

vid= VideoReader(fullFileName); %read the video
vid.NumberOfFrames

rgbImage1 = read(vid,1);
rgbImage = imcrop(rgbImage1,crop_drop);
grayImg_main = rgb2gray(rgbImage);

height = zeros(floor((frame_end-frame_begin)/frame_skip),size(grayImg_main,2));
time   = zeros(floor((frame_end-frame_begin)/frame_skip),1);

%starting looping over frames for analysis
ii = 1;
 for frame = frame_begin:frame_skip:frame_end %vid.NumberOfFrames;
all = [0 0 0 0]; %initializing the matrix with all the data
rgbImage1 = read(vid,frame);

%find the drop location now
rgbImage = imcrop(rgbImage1,crop_drop);
rgbImage_col=rgbImage;
bwImage=rgb2gray(rgbImage);

%figure()
%imshow(rgbImage)

%{
%% Step 1 Get Undistorted image
%Removing the distortion due to curvature effects of the camera
%Needs calibration
% URL: https://www.mathworks.com/help/vision/ug/remove-distortion-from-an-image-using-the-cameraparameters-object.html
% Better URL: https://conservancy.umn.edu/bitstream/handle/11299/173817/cameraParameters%20class.pdf?sequence=12&isAllowed=y
IntrinsicMatrix = [715.2699 0 0; 0 711.5281 0; 565.6995 355.3466 1];
radialDistortion = [-0.3361 0.0921]; 
cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion); 
J = undistortImage(rgbImage,cameraParams);
%figure; imshowpair(imresize(rgbImage,0.5),imresize(J,0.5),'montage');
%title('Original Image (left) vs. Corrected Image (right)');
%}

%% Step #2 Threshold
grayImage=rgb2gray(rgbImage);
grayImage = 255 - grayImage;
bwImage = im2bw(grayImage,1-70/255);

%figure()
%imshow(bwImage)



%% Step #3 Finding the surface
n=1;
k = ones(size(bwImage,2),1);
Ii = ones(size(bwImage,2),1);
Ib = ones(size(bwImage,2),1);
while n < size(bwImage,2)+1
    I = find(bwImage(:,n));
    %I = sort(I,'descend');
    if (size(I,1) > 10)
        for i=1:size(I,1)-5    %Checking the thickness of the white region to discard the grids 
            if (I(i)+1==I(i+1)&&I(i)+2==I(i+2)&&I(i)+3==I(i+3)&&I(i)+4==I(i+4)&&I(i)+5==I(i+5)&&I(i)+6==I(i+6))
                Ii(n) = I(i);
                break;
            else 
                Ii(n)=1;
            end
  
        end
        
        for i=size(I,1):-1:5    %Checking the thickness of the white region to discard the grids 
            if (I(i)-1==I(i-1)&&I(i)-2==I(i-2)&&I(i)-3==I(i-3)&&I(i)-4==I(i-4)&&I(i)-5==I(i-5)&&I(i)-6==I(i-6))
                Ib(n) = I(i);
                break;
            else 
                Ib(n)=1;
            end
        end
        
        
    end
    Ii(Ii<3)=1;
    %Ib(Ib<size(I,1)-100)=1;
    %k(n) = max(I);
    
    %if k(n)<5; k(n)=nan; end; %Naive approach to remove grid 
    rgbImage_col(Ii(n), n,:) = [255,0,0]; % or [255,255,255] if that doesn't work.
    rgbImage_col(Ib(n), n,:) = [0,255,0];
    n = n + 1;
    
end

height(ii,:) = Ib - Ii;
time(ii) = (frame-frame_begin)/fps;

%// Step #3 Find regions of drops
%rp = regionprops(im_thresh, 'BoundingBox', 'Area');

figure('Visible','Off')
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
subplot(1,1,1)
text_str = ['Time: ' num2str((frame-frame_begin)/fps,'%0.2f') 's'];
timer_img = insertText(rgbImage_col,position_timer,text_str,'FontSize',50,'BoxOpacity',0.0,'TextColor','white');
%imshow(rgbImage_col)
imshow(timer_img)
hold on
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
F(frame) = getframe(gcf) ;
hold off
drawnow


ii = ii+1;
 end

% create the video writer with fps of the original video
 Data_result= sprintf('%s_analyzed.mov',fffilename);
  writerObj = VideoWriter(Data_result);
  writerObj.FrameRate = fps; % set the seconds per image
  open(writerObj); % open the video writer
% write the frames to the video

for i=frame_begin:frame_skip:frame_end
    % convert the image to a frame
    frameimg = F(i) ;
    writeVideo(writerObj, frameimg);
end

% close the writer object
close(writerObj);
close all

outputfilename = append(fffilename,'_analysed','.mat');
save(outputfilename,'height','time');