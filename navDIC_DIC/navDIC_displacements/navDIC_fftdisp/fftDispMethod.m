function MovingPoints = fftDispMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)

    % PARALETERS
        CorrSize = 30 ;
        m = round(CorrSize/4) ; % Margin to truncate borders
        uMax = CorrSize/4 ; % Maximum allowed displacement per iteration
    
    % INFOS
        nPts = size(PtsMov,1) ;
        startTime = tic ;
        
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
        
    % Normalize
        fftRef = fftRef./abs(fftRef) ;
        fftMov = fftMov./abs(fftMov) ;
        
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
                [U,~,V] = svd(phi) ;
            % Shift Invariance
                k(p,2) = U(1:end-1,1)\U(2:end,1) ;
                k(p,1) = conj(V(1:end-1,1)\V(2:end,1)) ;
        end
        
    % Displacement
        U = real(1i*log(k))*CorrSize/2/pi ;
        normU = sum(U.^2,2) ;
        U(normU>uMax^2,:) = NaN ;
        
    % Moving Points
        MovingPoints = U + disPtsMov - (disPtsRef-PtsRef) ;
        
        
    % Timing
        disp(['fftdisp: ',num2str(toc(startTime)*1000,'%.1f'),' ms']) ;
end


% Dx = real(1i*log((wX.*UpX)\(wX.*DwnX)))*c/2/pi ;

        