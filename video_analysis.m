% Calculating the characteristics of a moving drop captured using a moving camera
% Code finds the number on the ruler using template matching (normxcorr2 algorithm), references the frame accordingly
% Code finds the drop and evaluates its characterstics considering parallax correction 
% Written by: Mohammad Afzal Shadab
% Last edited: May 11th, 2019
% Email: afzalshadab0@yahoo.com
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
threshold = 67; %setting the threshold of the image
crop_ruler= [1818 10 170 1700]; %cropping the ruler region: left, top, width, height
crop_drop=  [800 0 600 2000]; %cropping the drop region: left, top, width, height
frame_begin = 1;  %starting frame
frame_end = 78;  %end frame
position_timer = [0, 50];   %position of the timer
ii_start_thresh = 1;        %starting scale reading threshold
ii_end_thresh = 56;         %end scale reading threshold
ii_start=   1;              %first loop beginning over the template images
ii_end=     10;             %first loop ending over the template images
ii_thresh = 3;              %check for ii-ii_threshold to ii+ii_threhold templates (for optimization)

%Parallax referencing
S1=0.02301496;       %Scale at normal/ruler distance (mm/pixel)
S2=0.026954178;      %Scale at drop distance (mm/pixel)
parallax_ref=1000.0; %Reference pixel for no parallax
%freesurfpixel = 362;%Manually set the free surface pixel

%Read template files for matching
for ii = ii_start_thresh:ii_end_thresh
   ffolder = fullfile('Template');
   ffilename = sprintf('%d.png',ii);
   fullonFileName = fullfile(ffolder, ffilename);
   TemplateImage(:,:,:,ii) = imread(fullonFileName);
end

% Read the video in a standard MATLAB color video format
folder = fullfile('\');
fffilename = 'Sample';
baseFileName = sprintf('%s.avi',fffilename);
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

%starting looping over frames for analysis
 for frame = frame_begin:frame_end %vid.NumberOfFrames;
all = [0 0 0 0]; %initializing the matrix with all the data
rgbImage1 = read(vid,frame);
rgbImage = imcrop(rgbImage1, crop_ruler);
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(rgbImage);
for ii = ii_start:ii_end

  smallSubImage= TemplateImage(:,:,:,ii);

%smallSubImage = imread('D:\Research\MIT Research\Rheometer data\Image\1.png');
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(smallSubImage);

% Ask user which channel (red, green, or blue) to search for a match.
% channelToCorrelate = menu('Correlate which color channel?', 'Red', 'Green', 'Blue');
% It actually finds the same location no matter what channel you pick, 
% for this image anyway, so let's just go with red (channel #1).
% Note: If you want, you can get the template from every color channel and search for it in every color channel,
% then take the average of the found locations to get the overall best location.
channelToCorrelate = 1;  % Use the red channel.
correlationOutput = normxcorr2(smallSubImage(:,:,1), rgbImage(:,:, channelToCorrelate));
%subplot(2, 2, 3);
%imshow(correlationOutput, []);
%axis on;
% Get the dimensions of the image.  numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(correlationOutput);

% Find out where the normalized cross correlation image is brightest.
[maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
[yPeak, xPeak] = ind2sub(size(correlationOutput),maxIndex(1));
% Because cross correlation increases the size of the image, 
% we need to shift back to find out where it would be in the original image.
corr_offset = [(xPeak-size(smallSubImage,2)) (yPeak-size(smallSubImage,1))];

if (maxCorrValue>all(2) && corr_offset(2)>0) %Highest correlation template found in the image
    all = [ii maxCorrValue corr_offset(1) corr_offset(2) 100 1];
end
end
ii = all(1); 

%taking full image
crop_rulernew=crop_ruler;
crop_rulernew(4)=crop_rulernew(4)+300;
rgbImage = imcrop(rgbImage1, crop_rulernew);

%find the reference pixel on the scale
k=1000.0;
for i=all(4)+120:all(4)+200 %loop below the image
if rgbImage(i,all(3)+1,1)<k %all(3) is the x coordinate of bounding box
     y_ref = i;
     k = rgbImage(i,all(3),1);
end
end

%find the drop location now
rgbImage = imcrop(rgbImage1,crop_drop);
rgbImage=rgb2gray(rgbImage);
%// Step #1 Threshold
im_thresh = rgbImage < threshold;
interface = 0;
%// Find the interface if the scale value is less than 8cm
if ii<8;
    for j = 1:2000;
        if mean(im_thresh(j,:))>=0.4;
            im_thresh(1:j,:)=0;
            interface =1;
        end
    end
end

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