% Calculating the characteristics of a groundwater flow using a fixed camera
% Written by: Mohammad Afzal Shadab
% Last edited: June 20th, 2021
% Email: mashadab@utexas.edu
% Input: raw snapshot 
% Output: analysed image, datasheet
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
threshold = 40; %setting the threshold of the image

%Hydraulic conductivity measurement
W  = 2.54/100;%Width of the acrylic cell (m)
K  = 0.091;   %For 0.5mm= 0.0024 m/s; 1mm=0.091 m/s; 2mm=0.0285 m/s
h2 = 0;       %lake level (m)
%Manually setting the region of importance
top_height = 850;
bottom_height = 1125;
x_origin = 336;
x_right  = 4680;
scale_left = 322;
scale_right= 4586;
scale_length= 160; %length of the 16 boxes

L = (x_right-x_origin+1); %length in pixels

scale = scale_length/(scale_right-scale_left+1);  %conversion from pixels to height in cm (cm/pixel)
L = L*scale;
crop_drop=  [x_origin top_height (x_right - x_origin) (bottom_height - top_height)]; %cropping the drop region: left, top, width, height
frame_begin = 1;  %starting frame
frame_end   = 1;   %end frame
frame_skip = 60;

Q  = 25;     %Volumetric flow rate (mL/min)
Q  = Q*10^(-6)/60; %converting to m^3/s

% Read the video in a standard MATLAB color video format
folder = fullfile('\Images\1mm\');
fffilename = '25mL_per_minute_1mm_beads_8_June';
baseFileName = sprintf('%s.jpg',fffilename);
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

rgbImage1 = imread(fullFileName);
rgbImage = imcrop(rgbImage1,crop_drop);
grayImg_main = rgb2gray(rgbImage);

height = zeros(size(grayImg_main,2));

%starting looping over frames for analysis
ii = 1;
 for frame = frame_begin:frame_skip:frame_end %vid.NumberOfFrames;
all = [0 0 0 0]; %initializing the matrix with all the data
rgbImage1 = imread(fullFileName);

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
grayImage(grayImage<20)=255;
grayImage = 255 - grayImage;
bwImage = im2bw(grayImage,1-threshold/255);

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
    if (size(I,1) > 6)
        for i=1:size(I,1)-5    %Checking the thickness of the white region to discard the grids 
            if (I(i)+1==I(i+1)&&I(i)+2==I(i+2)&&I(i)+3==I(i+3)&&I(i)+4==I(i+4)&&I(i)+5==I(i+5))%&&I(i)+6==I(i+6)&&I(i)+7==I(i+7)&&I(i)+8==I(i+8)&&I(i)+9==I(i+9)&&I(i)+10==I(i+10)%&&I(i)+11==I(i+11))%&&I(i)+12==I(i+12)&&I(i)+13==I(i+13)&&I(i)+14==I(i+14)&&I(i)+15==I(i+15)&&I(i)+16==I(i+16))
                Ii(n) = I(i);
                break;
            else 
                Ii(n)=1;
            end
  
        end
        
        for i=size(I,1):-1:5    %Checking the thickness of the white region to discard the grids 
            if (I(i)-1==I(i-1)&&I(i)-2==I(i-2)&&I(i)-3==I(i-3))%&&I(i)-4==I(i-4)&&I(i)-5==I(i-5))
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

x     = linspace(1,size(height,2),size(height,2));
x     = (x-1)*scale;
height= (Ib - Ii+1)*scale;

TF = find(isoutlier(height,'movmedian',100));
height(TF) = NaN;

%// Step #3 Find regions of drops
%rp = regionprops(im_thresh, 'BoundingBox', 'Area');

figure('Visible','On')
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
subplot(1,1,1)
imshow(rgbImage_col)
hold on
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
F(frame) = getframe(gcf) ;
hold off
drawnow


ii = ii+1;
 end
xexp = x/100;
hexp = height/100;

outputfilename = append(fffilename,'_analysed','.mat');
save(outputfilename,'hexp','xexp','Q','W','L','K');

figure()
plot(xexp,hexp,'r.')
ylim([0,0.4])
xlim([0,1.7])
xlabel('x (m)')
ylabel('height (m)')
outputfilename = append(fffilename,'_analysed','.pdf');
saveas(gcf,outputfilename)