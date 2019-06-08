%% IF BACKWARD MODE, SOME POST-PROCESSING ARE NEEDED
    if reverseReference && dicDir<0
        % Backup Old Nodes
            oldNodes = Nodes ;
        % Reverse the displacements
            % First Valid Time
                firstValidFrame = sum(isnan(Xn(:,1,:)),3)+1 ;
            % New Nodes as corresponding position
                Nodes = [...
                        reshape(Xn(sub2ind(size(Xn),1:nNodes,ones(1,nNodes),firstValidFrame(:)')),[nNodes 1]) ...
                        reshape(Xn(sub2ind(size(Xn),1:nNodes,2*ones(1,nNodes),firstValidFrame(:)')),[nNodes 1]) ...
                                    ] ;
            % New Displacements
                Un = Xn-repmat(Nodes,[1 1 nFrames]) ;
        % Re-process the reference geometry with the new Nodes
            globalDIC_02_1_ProcessGeometries ;
    end

%% COMPUTE STRAINS
    % Valid elements as function of time
        VALID_t = reshape(~isnan(Un(:,1,:)),[nNodes nFrames]) ;
        validElems_t = false(nElems,nFrames) ;
        validElems_t(:,1:refFrame) = true ;
        for ii = 1:nFrames % no need to go to nImages if it has been stopped before
            validElems_t(:,ii) = sum(tri2nod(VALID_t(:,ii),:),1)'==3 ;
        end
    % Strains
        % DOFs as a function of time
            An = reshape(Un,[2*nNodes nFrames]) ; 
            An(isnan(An)) = 1i ; % Propagates NaNs with complex numbers
        % Gradient (at the element level)
            dUx_dx = gradT{1,1}*An + 1i*~validElems_t ;
            dUx_dy = gradT{1,2}*An + 1i*~validElems_t ;
            dUy_dx = gradT{2,1}*An + 1i*~validElems_t ;
            dUy_dy = gradT{2,2}*An + 1i*~validElems_t ;
        % Strains with NL terms
            Strains = [] ;
            Strains(:,1,:) = dUx_dx + 0.5*(dUx_dx.^2 + dUy_dx.^2) ;
            Strains(:,2,:) = dUy_dy + 0.5*(dUx_dy.^2 + dUy_dy.^2) ;
            Strains(:,3,:) = 0.5*(dUx_dy + dUy_dx + dUx_dx.*dUx_dy + dUy_dx.*dUy_dy) ;
            Strains(imag(Strains)~=0) = NaN ; % Re-set imaginary results to NaN
        % To the nodes level
            %meanOnElems = diag(1./sum(tri2nod(:,validElems),2))*tri2nod(:,validElems) ;
            
            
%% IF BACKWARD MODE, RESET TO THE ORIGINAL SEED MESH
    if reverseReference && dicDir<0
        % Retrieve Old Nodes
            Nodes = oldNodes ;
        % Re-process the reference geometry
            globalDIC_02_1_ProcessGeometries ;
    end
            
            
