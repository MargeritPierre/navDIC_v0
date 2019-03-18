function PtsApp = camsTo3d(hd,Pts)

Cam(1) = hd.Cameras(1).Properties ; 
Cam(2) = hd.Cameras(2).Properties ; 
PtsApp = zeros(size(Pts,1), 3, 2) ;

for j = 1:length(Cam)
    for i = 1 : size(Pts,1)
        PtsApp(:,:,j) = [ Cam(j).do / Cam(j).fx * ( Pts(i,1,j) - Cam(j).px / 2 ) ;...
                    Cam(j).do / Cam(j).fy * ( Pts(i,2,j) - Cam(j).py / 2 ) ; Cam(j).do ] ;
    end
end

end
