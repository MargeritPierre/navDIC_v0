function MovingPoints = fftDispMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)

    % INFOS
        nPts = size(PtsMov,1) ;
        
    % Discrete positions
        disPtsRef  = floor(PtsRef) ;
        disPtsMov  = floor(PtsMov) ;

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
        EXP = fftshift(fftshift(fftMov./fftRef,1),2) ;
    
    % Shift invariance
        upI = reshape(EXP(1:end-1,:,:),[],nPts) ;
        dwnI = reshape(EXP(2:end,:,:),[],nPts) ;
        upJ = reshape(EXP(:,1:end-1,:),[],nPts) ;
        dwnJ = reshape(EXP(:,2:end,:),[],nPts) ;
        kI = upI\dwnI ;
        kJ = upJ\dwnJ ;
end

        