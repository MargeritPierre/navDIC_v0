%% RUN GLOBAL DIC !

lastPlotTime = tic ;
for ii = dicFrames
    % Init
        it = 0 ;
        outFlag = false ;
    % Import and display current image
        img2 = Smooth(IMG(:,:,:,ii)) ;
        %im.CData = repmat(img2,[1 1 3]) ;
        im.CData = repmat((img2-min(img2(:)))/range(img2(:)),[1 1 3]) ;
    % First Guess for the displacement if needed
        if ~startWithNavDICPositions || all(isnan(Xn(:,1,ii)))
            % Take the previously computed frame
                Xn(:,:,ii) = Xn(:,:,ii-dicDir) ;
            % Add the previous "speed" as convergence help
                if addPreviousVelocity && abs(refFrame-ii)>=2
                    Xn(:,:,ii) = Xn(:,:,ii) + (Xn(:,:,ii-dicDir)-Xn(:,:,ii-2*dicDir)) * (frames(ii)-frames(ii-dicDir))/(frames(ii-dicDir)-frames(ii-2*dicDir)) * 1.0 ;
                end
            % Displacement
                Un(:,:,ii) = Xn(:,:,ii) - Nodes ;
        end
    % Newton-Raphson
        RMSE_0 = Inf ;
        while ~outFlag && ~stopBtn.Value
            % VALID GEOMETRY
                % Elements
                    validElems = validElems & sum(tri2nod(VALID,:),1)'==3 ; % triangles with still their three nodes valid
                    %validElems = validElems & sum(INSIDE(~deadPixels,:),1)'~=0 ; % triangles with still some pixels in the image
                % Edges
                    validEdges = validEdges & sum(edg2nod(VALID,:),1)'==2 ; % edges with still their two ends points valid
                    nakedEdges = sum(tri2edg(:,validElems),2)<2 ; % edges with less than two elements valid
                % Nodes
                    VALID = VALID & sum(tri2nod(:,validElems),2)>0 ; % nodes must be linked to at least one valid element
                    VALID = VALID & sum(edg2nod(:,validEdges),2)>0 ; % nodes must be linked to at least one valid edge
                    nVALID = sum(VALID) ;
                    if nVALID==0 ; break ; end
            % PIXEL-WISE DISPLACEMENT AND DIC DOMAIN
                % VALID DIC Domain
                    dicDomain = find(any(INSIDE(:,validElems),2)) ; 
                    selectDomain = sparse(1:length(dicDomain),dicDomain,1,length(dicDomain),nI*nJ) ;
                    dicMAPPING = selectDomain*MAPPING(:,VALID) ;
                % Compute the displacement 
                    Up = dicMAPPING*Un(VALID,:,ii) ;
                % New position of each pixel
                    JJp = JJ(dicDomain)+Up(:,1) ;
                    IIp = II(dicDomain)+Up(:,2) ;
                % Cull out-of-frame pixels
                    if ~cullOutOfFrame % pixels can be outside only if Nodes are outside
                        outOfFrame = IIp<1 | IIp>nI | JJp<1 | JJp>nJ ;
                        dicDomain(outOfFrame) = [] ;
                        JJp(outOfFrame) = [] ;
                        IIp(outOfFrame) = [] ;
                        selectDomain(outOfFrame,:) = [] ;
                        dicMAPPING(outOfFrame,:) = [] ;
                    end
            % IMAGE WARPING
                % Reduce the interpolation frame
                    iii = floor(min(IIp)):ceil(max(IIp)) ;
                    jjj = floor(min(JJp)):ceil(max(JJp)) ;
                % Warp the image
                    img2v = interp2(JJ(iii,jjj),II(iii,jjj),img2(iii,jjj),JJp,IIp,imWarpInterpOrder) ;
            % IMAGE MOMENTS
                % Select the part of the reference image of interest
                    img1v = img1(dicDomain) ;
                % Integration weights
                    switch size(localWEIGHT,2)
                        case nNodes
                            WEIGHT = localWEIGHT(:,VALID) ;
                        case nElems
                            WEIGHT = localWEIGHT(:,validElems) ;
                    end
                    WEIGHT = selectDomain*WEIGHT ;
                    sumWEIGHT = sum(WEIGHT,1).' ;
                % Mean over elements
                    meanImg1 = (WEIGHT'*img1v)./sumWEIGHT ;
                    meanImg2 = (WEIGHT'*img2v)./sumWEIGHT ;
                % Zero-local-mean on pixels
                    img1m = img1v-WEIGHT*meanImg1 ;
                    img2m = img2v-WEIGHT*meanImg2 ;
                % Norm over element
                    normImg1 = sqrt(WEIGHT'*(img1m.^2)) ;
                    normImg2 = sqrt(WEIGHT'*(img2m.^2)) ;
                % Zero-local-mean-normalized images
                    img1mz = img1m./(WEIGHT*normImg1) ;
                    img2mz = img2m./(WEIGHT*normImg2) ;
            % IMAGE FUNCTIONAL
                switch diffCriterion
                    case 'Diff' % Simple difference
                        diffImg = img1v-img2v ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_Diff' % Zero-mean difference
                        diffImg = img1m-img2m ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_N_Diff' % Normalized Zero-mean difference
                        diffImg = img1mz-img2mz ;
                        switch size(localWEIGHT,2)
                            case nNodes
                                ww = 1./normImg1(VALID) ;
                            case nElems
                                ww = 1./bsxfun(@(x,y)x./y,tri2nod(VALID,validElems)*normImg1,sum(tri2nod(VALID,validElems),2)) ;
                        end
                        weight = sparse(1:2*nVALID,1:2*nVALID,[ww;ww]) ;
                end
            % NEWTON-RAPHSON PROCEDURE
                % Recompute Jacobian and Hessian if needed
                    if refImageChanged
                        % Image gradients
                            dImg1_dx = dI_dx(img1) ;
                            dImg1_dy = dI_dy(img1) ;
                        % Jacobian
                            dImg1_da = [...
                                    spdiag(dImg1_dx(dicDomain))*dicMAPPING ...
                                    spdiag(dImg1_dy(dicDomain))*dicMAPPING ...
                                    ] ;
                        % Hessian
                            Hess = dImg1_da'*dImg1_da ;
                    end
                % Compute the first RMSE derivative
                    dr_da = (dImg1_da'*diffImg) ;
                % Contraint on the SECOND displacement gradient
                    validDOF = [VALID;VALID] ;
                    switch strainCriterion
                        case 'normal'
                            vEdg = repmat(~nakedEdges,[2 1]) ;
                            vEle = [validElems;validElems;validElems] ;
                            CONS = B(vEle,validDOF)'*(E(vEdg,vEle)'*E(vEdg,vEle))*B(vEle,validDOF) ;
                        case 'full'
                            vEdg = repmat(~nakedEdges,[4 1]) ;
                            vEle = [validElems;validElems;validElems;validElems] ;
                            CONS = G(vEle,validDOF)'*(E(vEdg,vEle)'*E(vEdg,vEle))*G(vEle,validDOF) ;
                    end
                    % Debugging...
                        %markers.XData = Xn(nodesOnNaked,1,ii); markers.YData = Xn(nodesOnNaked,2,ii);
                        %if any(~VALID)  disp('STOP!'); pause ; end
                % Updating DOFs, X*a=b
                    X = ( ...
                            weight*Hess*weight...
                            + beta*CONS...
                        ) ; 
                    b = (...
                            weight*dr_da...
                            - beta*CONS*[Un(VALID,1,ii);Un(VALID,2,ii)]...
                        ) ;
                    a = X\b ;
                % Displacement
                    Un(VALID,1,ii) = Un(VALID,1,ii) + a(1:nVALID) ;
                    Un(VALID,2,ii) = Un(VALID,2,ii) + a(nVALID+(1:sum(VALID))) ;
                % Positions
                    Xn(:,:,ii) = Nodes + Un(:,:,ii) ;
            % CONVERGENCE CITERIONS
                % Residues Maps
                    resid = abs(diffImg) ;
                    meanSquaredElemResid = sqrt((WEIGHT'*abs(diffImg).^2)) ;
                    corrCoeff = abs(WEIGHT'*(img1mz.*img2mz)) ;
                % Criterions
                    it = it+1 ;
                    %RMSE = norm(diffImg(DOMAIN))/norm(img1(DOMAIN)) ; 
                    normA = max(abs(a)) ; norm(a)/nNodes ;
                % Convergence criterion
                    if it>maxIt ; outFlag = true ; end
                    if normA<minNorm ; outFlag = true ; end
                    %if RMSE<1e-6 || abs((RMSE-RMSE_0)/RMSE) < 1e-4 ; outFlag = true ; end
                % Keep the error
                    %RMSE_0 = RMSE ;
            % POINTS/ELEMENTS VALIDATION/DELETION
                % OUT-OF-FRAME POINTS
                    if cullOutOfFrame
                        VALID = VALID & Xn(:,1,ii)<nJ+1 & Xn(:,1,ii)>0 & Xn(:,2,ii)<nI+1 & Xn(:,2,ii)>0 ;
                    end
                % DECORRELATED ELEMENTS
                    if (outFlag || alwaysValidGeometry)
                        switch size(localWEIGHT,2) 
                            case nElems % Correlation at the element level
                                %VALID = VALID & (tri2nod*corrCoeff(:))./sum(tri2nod(:,validElems),2) ;
                                cullElems = corrCoeff<minCorrCoeff | meanSquaredElemResid>maxMeanElemResidue ;
                                validElems(validElems) = validElems(validElems) & ~cullElems ; 
                            case nNodes % Correlation at the node level
                                VALID(VALID) = VALID(VALID) & corrCoeff(:)>minCorrCoeff ;
                        end
                    end
                % SET THE NON-VALID NODES TO NAN
                    Xn(~VALID,:,ii) = NaN ;
                    Un(~VALID,:,ii) = NaN ;
            % DISPLAY
                % Always display infos
                    infosText.String = [ ' INFOS' ...
                                         ' | Frame: ' num2str(frames(ii)) ...
                                         ' | It: ' num2str(it) ...
                                         ' | normA: ' num2str(normA,2) ...
                                         ] ;
                % Heavy Plots
                    if plotEachIteration || (outFlag && plotEachFrame) || toc(lastPlotTime)>1/plotRate % (outFlag && toc(lastPlotTime)>1/plotRate)
                        %ttl.String = [num2str(frames(ii)),'(',num2str(it),')'] ;
                        residues = ... WEIGHT*abs(WEIGHT'*(img1mz.*img2mz)) ... Correlation coeffient
                                   ... abs(diffImg) ... NR residues
                                    WEIGHT*sqrt((WEIGHT'*abs(diffImg).^2)) ... Sum of squared NR residues
                                   ;
                        imRes.CData = reshape(sparse(dicDomain,1,residues,nI*nJ,1),[nI nJ]) ;
                        mesh.Vertices = Xn(:,:,ii) ;
                        mesh.Faces = Elems(validElems,:) ;
                        %figure(figDebug) ; clf ; ind = (0:nJ-1)*nI+ceil(nI/2) ; plot(img1v(ind)) ; plot(img2v(ind)) ; 
                        lastPlotTime = tic ;
                        drawnow ;
                    end
            % Pause execution ?
                if pauseAtPlot && ~continueBtn.Value
                    while ~stopBtn.Value && ~nextBtn.Value && ~continueBtn.Value
                        drawnow ;
                    end
                    nextBtn.Value = 0 ;
                    if continueBtn.Value
                        nextBtn.Visible = 'off' ;
                        continueBtn.Visible = 'off' ;
                    end
                end
        end
        % Modify the reference image if needed
            if weightCurrentImage>0
                img1(dicDomain) = img1(dicDomain)*(1-weightCurrentImage) + img2v*weightCurrentImage ;
                refImageChanged = true ;
            end
        % Out Criterions
            if stopBtn.Value; break ; end
            if nVALID==0; break ; end
end