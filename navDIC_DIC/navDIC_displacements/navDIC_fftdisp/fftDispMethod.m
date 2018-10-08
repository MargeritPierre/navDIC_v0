function MovingPoints = fftDispMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)

    % INFOS
        CorrSize = 20 ;
        m = round(CorrSize/4) ; % Margin to truncate borders
        r = 1 ; % order of the FFT ;
        nPts = size(PtsMov,1) ;
        
    % Discrete positions
        disPtsRef  = round(PtsRef) ;
        disPtsMov  = round(PtsMov) ;

    % Create Imagettes
        % Indices
            % One imagette centered in 0
                ii = -floor(CorrSize/2)+(1:CorrSize) ;
                [JJ,II] = meshgrid(ii,ii) ;
                III = repmat(II,[1 1 nPts]) ;
                JJJ = repmat(JJ,[1 1 nPts]) ;
            % Imagette for each point
                IIIRef = bsxfun(@plus,III,reshape(disPtsRef(:,2),[1 1 nPts])) ;
                JJJRef = bsxfun(@plus,JJJ,reshape(disPtsRef(:,1),[1 1 nPts])) ;
                IIIMov = bsxfun(@plus,III,reshape(disPtsMov(:,2),[1 1 nPts])) ;
                JJJMov = bsxfun(@plus,JJJ,reshape(disPtsMov(:,1),[1 1 nPts])) ;
            % Linear Indices
                NNNRef = sub2ind(size(imgRef),IIIRef,JJJRef) ;
                NNNMov = sub2ind(size(imgMov),IIIMov,JJJMov) ;
            % Pixel Values
                imagettesRef = imgRef(NNNRef) ;
                imagettesMov = imgMov(NNNMov) ;
    
    % FFT
        fftRef = fft(fft(imagettesRef,[],1),[],2) ;
        fftMov = fft(fft(imagettesMov,[],1),[],2) ;
        
    % Phase field
        PHI = fftshift(fftshift(fftMov./fftRef,1),2) ;
        PHI = PHI./abs(PHI) ;
        PHI(isnan(PHI)) = 1 ;
    
    % Shift invariance
%         upI = reshape(PHI(1:end-1,:,:),[],nPts) ;
%         dwnI = reshape(PHI(2:end,:,:),[],nPts) ;
%         upJ = reshape(PHI(:,1:end-1,:),[],nPts) ;
%         dwnJ = reshape(PHI(:,2:end,:),[],nPts) ;
        
    % Wavevectors
        k = zeros(nPts,2) ;
        for p = 1:nPts
            % Select and truncate
                phi = PHI(1+m:end-m,1+m:end-m,p) ;
            % SVD Filtering
                [U,S,V] = svd(phi) ;
                phi = U(:,1:r)*S(1:r,1:r)*V(:,1:r)' ;
            % Shift Invariance
                upI = reshape(phi(1:end-1,:),[],1) ;
                dwnI = reshape(phi(2:end,:),[],1) ;
                upJ = reshape(phi(:,1:end-1),[],1) ;
                dwnJ = reshape(phi(:,2:end),[],1) ;
            % Wavevector estimation
                k(p,2) = upI\dwnI ;
                k(p,1) = upJ\dwnJ ;
        end
        
    % Displacement
        U = real(1i*log(k))*CorrSize/2/pi ;
        MovingPoints = U + disPtsMov ;
        
end


% Dx = real(1i*log((wX.*UpX)\(wX.*DwnX)))*c/2/pi ;

        