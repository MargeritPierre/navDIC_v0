OT_Corrected_Data = hd.AcquiredOTPositions.Data;
for i=1:size(hd.AcquiredOTPositions.Data,1)
    for k = 1:3
        for j=308:514
            OT_Corrected_Data(i,k,j) = hd.AcquiredOTPositions.Data(i,k,308) ...
                +(hd.AcquiredOTPositions.Data(i,k,514)-hd.AcquiredOTPositions.Data(i,k,308))...
                /(hd.AcquiredOTPositions.Time(514)-hd.AcquiredOTPositions.Time(308))...
                *(hd.AcquiredOTPositions.Time(j)-hd.AcquiredOTPositions.Time(308));
        end
    end
end