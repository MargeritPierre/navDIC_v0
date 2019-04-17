function PtsApp = camsTo3d(hd,Pts)
nbCam = length(hd.Cameras) ;
for i = 1:nbCam
    Cam(i) = hd.Cameras(i).Properties ; 
end
PtsApp = zeros(size(Pts,1), 3, nbCam) ;

for j = 1:length(Cam)
    for i = 1 : size(Pts,1)
        PtsApp(i,:,j) = [ Cam(j).do / Cam(j).fpix(1) * ( Pts(i,1,j) - Cam(j).PPL(1) ) ,...
                    Cam(j).do / Cam(j).fpix(2) * ( Pts(i,2,j) - Cam(j).PPL(2) ) , Cam(j).do ] ;
    end
end

end
