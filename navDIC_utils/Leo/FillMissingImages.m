% This scripts add up missing images into the NavDIC hd setup.
% I made this because of problems with data acquisition with a allied
% vision camera. Some frames were missing during acquisition, and when i'm
% doing the processing if the frame cell is empty, the fftdisp computation
% crashes.

% So the idea is to add the last valid image instead of the empty cell into
% the images cell array to avoid crashing of the computation. This is done
% using a simple for loop and duplicates the last valid image into the
% empty cells.

% Improvement to do : add a flag to identify duplicated frames to avoid
% misinterpretation of results.

global hd

if isempty(hd.Images{1, 1}{1, 1})
hd.Images{1, 1}{1, 1} = hd.Images{1, 1}{1, 3};
end

for i = 1:length(hd.Images{1, 1})
if isempty(hd.Images{1, 1}{1, i})
hd.Images{1, 1}{1, i} = hd.Images{1, 1}{1, i-1};
disp(strcat('Filled empty frame at index  ',num2str(i)));
end
end



