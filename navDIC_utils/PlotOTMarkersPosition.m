%THIS SCRIPT ENABLES TO PLOT THE EVOLUTION OF THE OT MAKERS POSITIONS
%THROUGHOUT THE ACQUISITION

%% DATA PARAMETERS

global hd;
OTPos = [];

%% CREATE MATRIX TO PLOT

time = floor(hd.AcquiredOTPositions.Time(1));
nbOfPositionsSelected = 0;
for i = 1:size(hd.AcquiredOTPositions.Time,1)
    if not(floor(hd.AcquiredOTPositions.Time(i)) == time)
        nbOfPositionsSelected = 0;
        time = floor(hd.AcquiredOTPositions.Time(i));
    end
    if nbOfPositionsSelected == 0
        OTPos = cat(3, OTPos,...
            hd.AcquiredOTPositions.Data(:,:,i));
        %ones(size(hd.AcquiredOTPositions.Data,1),1)* hd.AcquiredOTPositions.Time(i)
        nbOfPositionsSelected = nbOfPositionsSelected + 1;
    end
end

%% MEASURE DISPLACEMENT OVER TIME

d = [];
for i = 1:size(OTPos,3)
tempD = [];
    for j = 1:size(OTPos,1)
        tempD = [tempD pdist([OTPos(j,:,i) ; OTPos(j,:,1)],'euclidean')];
    end
    d = [d; tempD];
end
d = d - min(d);
d = d./max(d);

%% GENERATE PLOT

f1 = figure;
%set(gca, 'DataAspectRatio', [1,1,1]);
xlabel('x(m)')
ylabel('z(m)')
for i = 1:size(OTPos,1)
    c = [zeros(size(d,1),1) zeros(size(d,1),1) d(:,i)];
    scatter(OTPos(i,3,:),OTPos(i,2,:),50,c,'.')
end