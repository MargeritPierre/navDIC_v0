
global hd
clearvars -except hd

figure

%%

seed = 2 ;

camID = hd.Seeds(seed).CamIDs ;

IMG = hd.Images{camID}(:,:,:,hd.Seeds(seed).RefFrame) ;
%IMG = double(IMG) ;
%IMG = log10(IMG+1) ;

%IMG = medfilt2(IMG,20*[1 1]) ;

N = 5 ;
disk = sqrt((-N/2:N/2)'.^2 + (-N/2:N/2).^2)<N/2 ;
rect = true(N) ;
mask = rect ; N2 = sum(mask(:)) ;
%IMG = ordfilt2(IMG,N2,mask) ; % Maximum
%IMG = ordfilt2(IMG,1,mask) ; % Minimum
%IMG = ordfilt2(IMG,ceil(N2/2),mask) ; % Median
IMG = ordfilt2(IMG,N2,mask) - ordfilt2(IMG,1,mask)  ; % Range

filt = blackman(100) ;
IMG = conv2(IMG,filt(:)*filt(:)','same') ;

%gradX = conv2(IMG,[-1 0 1]/2,'valid') ;
%gradY = conv2(IMG,[-1 0 1]'/2,'valid') ;
%IMG = sqrt(gradX(2:end-1,:).^2 + gradY(:,2:end-1).^2) ;
%IMG = conv2(IMG,filt(:)*filt(:)','same') ;

IMG = double(IMG) ;
IMG = (IMG-min(IMG(:)))/range(IMG(:)) ;

clf
imagesc(IMG)
set(gca,'ydir','reverse')
axis equal
axis tight
axis off
colorbar
colormap gray




