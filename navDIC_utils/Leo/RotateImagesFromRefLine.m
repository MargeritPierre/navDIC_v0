% This script helps to get the X and Y of image well-aligned with some
% reference features in the image, such as a plane. This might be important because one
% can need to check the u1 (x) and u2 (y) displacements relatively to a
% plane.


%% Load reference image from hd
global hd

refImg = hd.Images{1, 1}{1, 1};

%% Use drawing tool to let user define the two points or draw line.

angle = DefineScale(refImg);

%% Load the two points

% Operate the angle

% A = [50  70]; %[x y] & x,y: pixel number
% B = [100 30]; %[x y]

A = point1;
B = point2;

% find point2 (point with larger x)
if A(1)>B(1)
  P2 = A; P1 = B;
else
  P2 = B; P1 = A;
end
% clockwise OR counterclockwise
if P2(2)<P1(2)
  w=-1; % clockwise
else
  w=+1; % counterclockwise
end
angle = w*atan(abs(P2(2)-P1(2))/(P2(1)-P1(1)))*180/pi; % deg

% Apply rotation to all images of hd.

ImgRot = imrotate(refImg,angle);
figure;
title("rotated")
imshow(ImgRot);

%% Go and rotate all images in navDIC memory.

for i = 1:length(hd.Images{1, 1})
hd.Images{1, 1}{1, i} = imrotate(hd.Images{1, 1}{1, i},angle);
disp(strcat('Rotated frame at index  ',num2str(i)));
end

%% Subfunctions

function angle = DefineScale(img)
figure;
a = imshow(img);
% start drawing
h = imdistline;
api = iptgetapi(h);
myvar = 0;
angle = api.getAngleFromHorizontal();
%el = addlistener(h,@)
fprintf('The angle of the segment is: %0.2f ??? \n', angle)
% while myvar == 0
% if myvar == 1
%     angle = api.getAngleFromHorizontal();
% %el = addlistener(h,@)
% fprintf('The angle of the segment is: %0.2f ??? \n', angle)
% end
% end
end