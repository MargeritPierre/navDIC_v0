function PtsApp = camsTo3d(hd,Pts)
nbCam = length(hd.Cameras) ;
for i = 1:nbCam
    Cam(i) = hd.Cameras(i).Properties ; 
end
PtsApp = zeros(size(Pts,1), 3, nbCam) ;

for j = 1:length(Cam)
    for i = 1 : size(Pts,1)
        PtsApp(i,:,j) = [ Cam(j).do / Cam(j).fx * ( Pts(i,1,j) - Cam(j).px/2 ) ,...
                    Cam(j).do / Cam(j).fy * ( Pts(i,2,j) - Cam(j).py/2 ) , Cam(j).do ] ;
    end
end

end
