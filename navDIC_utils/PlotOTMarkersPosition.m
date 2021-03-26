%THIS SCRIPT ENABLES TO PLOT THE EVOLUTION OF THE OT MAKERS POSITIONS
%THROUGHOUT THE ACQUISITION

%% DATA PARAMETERS

global hd2;
OTPos = [];

%% CREATE MATRIX TO PLOT

time = floor(hd2.AcquiredOTPositions.Time(1));
nbOfPositionsSelected = 0;
for i = 1:size(hd2.AcquiredOTPositions.Time,1)
    if not(floor(hd2.AcquiredOTPositions.Time(i)) == time)
        nbOfPositionsSelected = 0;
        time = floor(hd2.AcquiredOTPositions.Time(i));
    end
    if nbOfPositionsSelected == 0
        OTPos = cat(3, OTPos,...
            hd2.AcquiredOTPositions.Data(:,:,i));
        %ones(size(hd2.AcquiredOTPositions.Data,1),1)* hd2.AcquiredOTPositions.Time(i)
        nbOfPositionsSelected = nbOfPositionsSelected + 1;
    end
end

%% MEASURE DISPLACEMENT OVER TIME

d = [];
for i = 2:size(OTPos,3)
tempD = [];
    for j = 1:size(OTPos,1)
        %testMat = [OTPos(j,3,i) OTPos(j,2,i); OTPos(j,3,1) OTPos(j,2,1)]
        dist2D = pdist([OTPos(j,3,i) OTPos(j,2,i); OTPos(j,3,1) OTPos(j,2,1)]);
        %tempD = [tempD pdist([OTPos(j,:,i) ; OTPos(j,:,1)])];
        tempD = [tempD dist2D];
    end
    d = [d; tempD];
end
dNorm = d - min(d);
dNorm = dNorm./max(dNorm);

%% GENERATE PLOT

f1 = figure;
set(gca, 'DataAspectRatio', [1,1,1]);
xlabel('x(m)')
ylabel('z(m)')
title('Optitrack Marker Positions \& 2D Displacement');
subtitle('! Only 2D displacement are computed and plotted on this graph','FontSize',14,'Color',[0.87 0.34 0.28]);
lightBlue = [104 188 227]/255;
darkBlue = [11 82 115]/255;
for i = 1:size(OTPos,1)
    scatter(OTPos(i,3,1),OTPos(i,2,1),50,'red','*')
    c = [lightBlue(1) + (darkBlue(1) - lightBlue(1))* dNorm(:,i)...
        lightBlue(2) + (darkBlue(2) - lightBlue(2))* dNorm(:,i)...
        lightBlue(3) + (darkBlue(3) - lightBlue(3))* dNorm(:,i)];
    scatter(OTPos(i,3,2:size(OTPos,3)),OTPos(i,2,2:size(OTPos,3)),50,c,'.')
    [M, I] = max(d(:,i));
    text(double(OTPos(i,3,I+1)),double(OTPos(i,2,I+1)),strcat(num2str(d(I,i)*1000,'%.1f'),' mm $\rightarrow$'), 'HorizontalAlignment', 'right','FontSize',8);
end