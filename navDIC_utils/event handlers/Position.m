function Position( src , evnt )

freq = 1 ; % 
global lastOTTime
if isempty(lastOTTime) ; lastOTTime = tic ; end 
if toc(lastOTTime)<1/freq ; return ; end
lastOTTime = tic ;
disp('OTPositions acquired!')

    % Marker positions acquisition
    global hd;
    tempMatrixOfPositions = [];
    tempMatrixOfID = [];
%     tempMatrixOfPositions =    [[evnt.data.LabeledMarkers.x]' ...
%                                 [evnt.data.LabeledMarkers.y]' ...
%                                 [evnt.data.LabeledMarkers.z]'];
    for i = 1:hd.OTMarkersCount
        if ~isempty(evnt.data.LabeledMarkers(i))
            tempMatrixOfPositions = [tempMatrixOfPositions; ...
                                            [evnt.data.LabeledMarkers(i).x ...
                                             evnt.data.LabeledMarkers(i).y ...
                                             evnt.data.LabeledMarkers(i).z]];
            tempMatrixOfID = [tempMatrixOfID; ...
                                evnt.data.LabeledMarkers(i).ID];
        else
            tempMatrixOfPositions = [tempMatrixOfPositions; ...
                                            [0.0 0.0 0.0]];
            tempMatrixOfID = [tempMatrixOfID; ...
                                evnt.data.LabeledMarkers(i).ID];
        end
    end
    hd.AcquiredOTPositions.Data = cat(3,hd.AcquiredOTPositions.Data,...
                                            tempMatrixOfPositions);
    hd.AcquiredOTPositions.ID = cat(3,hd.AcquiredOTPositions.ID,...
                                            tempMatrixOfID);
    hd.AcquiredOTPositions.Time = [hd.AcquiredOTPositions.Time ; ...
                                            evnt.data.fTimestamp];
                                        
end