%% THIS SCRIPT GENERATES THE IMAGES OF A RIGID MOTION
figure

%% REFERENCE FRAME

    % Parameters
        sz = 512*[1 1] ; % Image size in pixels [nI nJ]
    	N = round(7/10*prod(sz)); % Number of initial points
        medFiltSize = 7 ; % median fliter size
        gaussFiltSize = 9 ; % gaussian fliter size
        frame = round([.2 .2 ; .6 .6] .* sz) ; % object containing rectangle [i0 j0 ; nI ; nJ]
        edgeWidth = 3 ;

    % Texture generation
        I0 = false(sz) ;
        % Random distribution of active pixels
            ii = randi([1 sz(1)],N,1) ;
            jj = randi([1 sz(2)],N,1) ;
            I0(sub2ind(size(I0),ii,jj)) = true ;
        % Median Filter
            I0 = medfilt2(I0,[1 1]*medFiltSize) ;
            
    % Object framing
        [JJ,II] = meshgrid(1:sz(2),1:sz(1)) ;
        inside = II>=frame(1,1) & II<=frame(1,1)+frame(2,1) & JJ>=frame(1,2) & JJ<=frame(1,2)+frame(2,2) ;
        I0(~inside) = 0 ;
        I0(inside & ~bwmorph(inside,'thin',edgeWidth)) = 1 ;
    
    % Gaussian filter
        f = gausswin(gaussFiltSize) ; f = f(:)*f(:)' ; f = f./sum(abs(f(:))) ;
        I0 = conv2(double(I0),f,'same') ;
        
    % Display
        clf
        imagesc(repmat(I0,[1 1 3]))
        axis tight
        axis equal
        axis ij
        box on
        
        
        
%% OBJECT MOTION

    % Parameters
        R = linspace(0,2*pi,361) ; % Rotation at each step
        T = [0 0] ; % Translation at each step (in pixels)
        
    % Generate steps
        % Format position vectors
            R = R(:) ;
            R = R + T(:,1)*0 ;
            T = T + [R R]*0 ;
        % Initialize
            nImages = length(R) ;
            I = zeros([sz nImages]) ;
        % Build Images
            for tt = 1:nImages
                I(:,:,tt) = imtranslate(imrotate(I0,R(tt).*180/pi,'bilinear','crop'),T(tt,:),'bilinear') ;
            end
        
    % Display
        clf
        im = imagesc(repmat(I0,[1 1 3])) ;
        axis tight
        axis equal
        axis ij
        box on
        for tt = 1:nImages
            im.CData = repmat(I(:,:,tt),[1 1 3]) ;
            drawnow ;
        end
        
        
%% SAVE IMAGES

    % Parameters
        commonName = 'img' ;

    [path] = uigetdir('Images','Select a folder where to save the images') ;
    if path==0 ; return ; end
        
    for tt = 1:nImages
        imwrite(I(:,:,tt),[path '/' commonName '_' num2str(tt-1) '.tif']) ;
    end
    
    disp('Images Saved!') ;





