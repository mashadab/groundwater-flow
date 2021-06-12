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
fps = 5; %frames rate (per second)
threshold = 70; %setting the threshold of the image

%fprintf("Hi");

%error('Breaking out of function1');


top_height = 321;
bottom_height = 543;
x_origin = 96;
x_right  = 1243;

crop_drop=  [x_origin top_height (x_right - x_origin) (bottom_height - top_height)]; %cropping the drop region: left, top, width, height
frame_begin = 4000;  %starting frame
frame_end = 4010;   %end frame
frame_skip = 1;
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
grayImg_main  = rgb2gray(rgbImage);

%starting looping over frames for analysis
 for frame = frame_begin:frame_begin%frame_skip:frame_end %vid.NumberOfFrames;
all = [0 0 0 0]; %initializing the matrix with all the data
rgbImage1 = read(vid,frame);

%find the drop location now
rgbImage = imcrop(rgbImage1,crop_drop);
rgbImage_col=rgbImage;
bwImage=rgb2gray(rgbImage);
imshow(rgbImage)
figure()
imshow(rgbImage)

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
%grayImage(grayImage < 30)= 255;  %Getting rid of grids
%bwImage(bwImage < 70)= 0;
%bwImage(bwImage >= 70)= 255;
grayImage = 255 - grayImage;
bwImage = im2bw(grayImage,1-70/255);
%im_thresh = rgbImage < threshold;
figure()
imshow(bwImage)



%% Step #3 Finding the surface
n=1;
k = ones(size(bwImage,2),1);
Ii = ones(size(bwImage,2),1);
while n < size(bwImage,2)+1
    I = find(bwImage(:,n));
    %I = sort(I,'descend');
    if (size(I,1) > 10)
        for i=1:size(I,1)-5    %Checking the thickness of the white region to discard the grids 
            if (I(i)+1==I(i+1)&&I(i)+2==I(i+2)&&I(i)+3==I(i+3)&&I(i)+4==I(i+4)&&I(i)+5==I(i+5)&&I(i)+6==I(i+6))
                Ii(n) = I(i)
                'working'
                break;
            else 
                Ii(n)=1;
            end
        end
    end
    Ii(Ii<3)=1;
    %k(n) = max(I);
    
    %if k(n)<5; k(n)=nan; end; %Naive approach to remove grid 
    rgbImage_col(Ii(n), n,:) = [255,0,0]; % or [255,255,255] if that doesn't work.
    n = n + 1;
end


%// Step #3 Find regions of drops
%rp = regionprops(im_thresh, 'BoundingBox', 'Area');

figure('Visible','On')
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
 end

%{
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

%}
 
 
 
%{
%// Step #3 Find regions of drops
rp = regionprops(im_thresh, 'BoundingBox', 'Area');

%if no drop found
if isempty(rp) %& frame ~= frame_end
    %frame = frame +1 
figure('Visible','Off')
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
subplot(2,2,[1 3])
text_str = ['Time: ' num2str((frame-frame_begin)/fps,'%0.2f') 's'];
timer_img = insertText(rgbImage1,position_timer,text_str,'FontSize',100,'BoxOpacity',0.0,'TextColor','black');
imshow(timer_img)
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
F(frame) = getframe(gcf) ;
hold on

subplot(2,2,2)
yyaxis left
loglog(data(:,1),data(:,8),'^')
xlim([0,(frame_end-frame_begin)/fps])
ylim([0,600])
ylabel('$$Height\quad(mm)$$','Interpreter','latex')
xlabel('$$Time\quad (s)$$','Interpreter','latex')
hold on
yyaxis right
loglog(data(:,1),data(:,9),'o')
 xlim([0,(frame_end-frame_begin)/fps])
%ylim([0,200])
 ylabel('$$Velocity\quad(mm/s)$$','Interpreter','latex')
 xlabel('$$Time\quad (s)$$','Interpreter','latex')
 hold on
F(frame) = getframe(gcf) ;

 subplot(2,2,4)
 semilogx(data(:,1),data(:,6),'^')
%semilogx(data(:,1),data(:,6)./data(:,7),'^')
xlim([0,(frame_end-frame_begin)/fps])
%ylim([0,3])
ylabel('$$Width \quad (mm)$$','Interpreter','latex')
xlabel('$$Time\quad (s)$$','Interpreter','latex')
hold on
F(frame) = getframe(gcf) ;
drawnow
hold off

    data= [data;(frame-frame_begin)/fps,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]; %Time (s), ruler reading (mm), reference pixel,bounding box - left pixel, top pixel, width, height 
    continue
end

%// Step #4 Calculating the parameters and plotting regions on the image
area = [rp.Area].';
[~,ind] = max(area);
bb = rp(ind).BoundingBox;
drop_end=bb(2)+bb(4);

if bb(4)>3.0*bb(3) %if height is greater than 3 times the width
    drop_top = drop_end-3.0*bb(3);
    bb(2) = drop_top; 
    bb(4)=3.0*bb(3);
end

area_drop=imcrop(im_thresh,bb);
area_drop=S2*S2*bwarea(area_drop); %area in mm^2

% Display the original color image.
boxRectruler = [all(3)+crop_ruler(1) all(4) 100 120];
% Plot the box over the image.    
figure('Visible','Off')
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);

subplot(2,2,[1 3])
text_str = ['Time: ' num2str((frame-frame_begin)/fps,'%0.2f') 's'];
timer_img = insertText(rgbImage1,position_timer,text_str,'FontSize',100,'BoxOpacity',0.0,'TextColor','black');
imshow(timer_img)
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
rectangle('position', boxRectruler, 'edgecolor', 'g', 'linewidth',2);
% Give a caption above the image.

%reference line
boxRectrefline = [all(3)+crop_ruler(1) y_ref 100 1]; 
rectangle('position', boxRectrefline, 'edgecolor', 'r', 'linewidth',1);

bb(1) = bb(1)+crop_drop(1); %correcting for the drop location
rectangle('Position', bb, 'EdgeColor', 'red');

%calculate actual height and velocity
height=ii*10.0+(bb(2)+bb(4)-boxRectrefline(2))*S1+(S2-S1)*(bb(2)+bb(4)-parallax_ref);
if frame == frame_begin
    velocity = 0.0; %setting the initial velocity to be zero
    freesurf=ii*10.0+(bb(2)-boxRectrefline(2))*S1+(S2-S1)*(bb(2)-parallax_ref); %free surface location
    %freesurf=ii*10.0+(freesurfpixel-boxRectrefline(2))*S1+(S2-S1)*(freesurfpixel-parallax_ref); %free surface location
elseif frame > frame_begin && frame < frame_begin+10
    velocity = (height-data(frame-frame_begin,8))/((frame-frame_begin)/fps); %for varying time frames
else
    velocity = (height-data(frame-frame_begin-9,8))/(10/fps); % for ten time steps
end
data= [data;(frame-frame_begin)/fps,ii*10.0,boxRectrefline(2), bb(1),bb(2),bb(3)*S2,bb(4)*S2,height,velocity,area_drop,height-freesurf,velocity/(bb(3)*S2)]; %Time (s), ruler reading (mm), reference pixel,bounding box - left pixel, top pixel, width, absolute height(mm), velocity (mmps),height from free surface (mm), shear rate = velocity/width (1/s) 

hold on

subplot(2,2,2)
yyaxis left
loglog(data(:,1),data(:,8),'^')
xlim([0,(frame_end-frame_begin)/fps])
ylim([0,600])
ylabel('$$Height\quad(mm)$$','Interpreter','latex')
xlabel('$$Time\quad (s)$$','Interpreter','latex')
hold on
yyaxis right
loglog(data(:,1),data(:,9),'o')
 xlim([0,(frame_end-frame_begin)/fps])
%ylim([0,200])
 ylabel('$$Velocity\quad(mm/s)$$','Interpreter','latex')
 xlabel('$$Time\quad (s)$$','Interpreter','latex')
 hold on
 
F(frame) = getframe(gcf) ;

 subplot(2,2,4)
 semilogx(data(:,1),data(:,6),'^')
%semilogx(data(:,1),data(:,6)./data(:,7),'^')
xlim([0,(frame_end-frame_begin)/fps])
%ylim([0,3])
ylabel('$$Width \quad (mm)$$','Interpreter','latex')
xlabel('$$Time\quad (s)$$','Interpreter','latex')

hold on
 
F(frame) = getframe(gcf) ;
drawnow
hold off

%resetting the template sample for matching in the next loop for
%optimization
 ii_start = all(1)-ii_thresh;
 ii_end = all(1)+ii_thresh;
 %setting the boundary of template sample for matching
 if ii_start<ii_start_thresh ii_start=ii_start_thresh; end  
 if ii_end>ii_end_thresh ii_end=ii_end_thresh; end 
 end
 
 % create the video writer with fps of the original video
 Data_result= sprintf('%s_analyzed.avi',fffilename);
  writerObj = VideoWriter(Data_result);
  writerObj.FrameRate = fps; % set the seconds per image
  open(writerObj); % open the video writer
% write the frames to the video
for i=frame_begin:frame_end
    % convert the image to a frame
    frameimg = F(i) ;
    writeVideo(writerObj, frameimg);
end
% close the writer object
close(writerObj);
close all

data_table = array2table(data,'VariableNames',{'Time_s', ...
        'Scaleref_mm','Scalerefpixel','boundingboxleft','boundingboxtop',...
        'Width','Height','Actualheight_mm','Velocity_mmps','Area_mm2','Heightfreesurf_mm','Shearrate_1pers'});
 
 Data_result= sprintf('%s.mat',fffilename);
 save(Data_result,'data_table')
%}
 