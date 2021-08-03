p = 0.015;
OTCorrected = zeros(size(hd.AcquiredOTPositions.Data,1),...
                    size(hd.AcquiredOTPositions.Data,2),...
                    size(hd.AcquiredOTPositions.Data,3));
OTCorrected(:,1,1) = hd.AcquiredOTPositions.Data(:,1,1);
OTCorrected(:,3,1) = hd.AcquiredOTPositions.Data(:,3,1);
OTCorrected(:,2,1) = hd.AcquiredOTPositions.Data(:,2,1);
OTCorrected(:,1,2) = hd.AcquiredOTPositions.Data(:,1,2);
OTCorrected(:,3,2) = hd.AcquiredOTPositions.Data(:,3,2);
OTCorrected(:,2,2) = hd.AcquiredOTPositions.Data(:,2,2);
% for i=1:size(hd.AcquiredOTPositions.Data,1)
%     for k = 1:3
%         m = 1;
%         referenceZ = hd.AcquiredOTPositions.Data(i,k,1);
%         for j=3:size(hd.AcquiredOTPositions.Data,3)
%             if(abs(hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,j-1))>0.0002) %%distance in mm
%                 referenceZ = OTCorrected(i,k,j-1);
%                 OTCorrected(i,k,j) = referenceZ; %+ p
%                 m = j;
%             else
%                 if(abs(hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,j-2))>0.0002)
%                     OTCorrected(i,k,j) = hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,m) + referenceZ;
%                 else
%                     OTCorrected(i,k,j) = hd.AcquiredOTPositions.Data(i,k,j);
%                 end  
%             end
%         end
%     end
% end
for i=1:size(hd.AcquiredOTPositions.Data,1)
    for k = 1:3
        m = 1;
        referenceZ = hd.AcquiredOTPositions.Data(i,k,1);
        for j=3:size(hd.AcquiredOTPositions.Data,3)
            if(abs(hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,j-1))>0.0002) %%distance in mm
                referenceZ = OTCorrected(i,k,j-1);
                OTCorrected(i,k,j) = referenceZ; %+ p
                m = j;
            else
                %if(abs(hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,j-2))>0.0002)
                    OTCorrected(i,k,j) = hd.AcquiredOTPositions.Data(i,k,j) - hd.AcquiredOTPositions.Data(i,k,m) + referenceZ;
                %else
                %    OTCorrected(i,k,j) = hd.AcquiredOTPositions.Data(i,k,j);
                %end  
            end
        end
    end
end