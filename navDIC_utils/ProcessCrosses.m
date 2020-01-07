%% PROCESS CROSSES TO COMPUTE MEMBERS ANGLES, ELONGATIONS, ETC.

% Parameters
    crossSeedNumber = 4 ;
    macroSeedNumber = 2 ;
    nMembers = 4 ;
    
    
% Retrieve data
    global hd
    Pos = hd.Seeds(crossSeedNumber).MovingPoints ;
    
% Reshape
    Pos = reshape(Pos,[],nMembers+1,2,hd.nFrames) ;
    Pos = permute(Pos,[1 3 4 2]) ;
    nNodes = size(Pos,1) ;
    
% To Edges
    Edg = Pos(:,:,:,1:nMembers)-Pos(:,:,:,nMembers+1) ;
% Complex notation
    Edg = Edg(:,1,:,:) + 1i*Edg(:,2,:,:) ;
    
% Edge Lengths
    EdgLen = abs(Edg) ;
    meanEdgLen = mean(EdgLen,4) ;
    devEdgLen = mean(abs(EdgLen-meanEdgLen),4) ;
% Edge Rotations
    EdgRot = angle(Edg(:,:,:,:)./Edg(:,:,1,:)) ;
    meanEdgRot = mean(EdgRot,4) ;
    devEdgRot = mean(abs(EdgRot-meanEdgRot),4) ;
    
% Compute cross seed fields, then add the new fields
    hd.Seeds(crossSeedNumber).computeDataFields ;
    hd.Seeds(crossSeedNumber).DataFields.CrossAnalysis = 'Crosses Analysis' ;
    hd.Seeds(crossSeedNumber).DataFields.MeanEdgLen = repmat(meanEdgLen,[nMembers+1 1 1]) ;
    hd.Seeds(crossSeedNumber).DataFields.DevEdgLen = repmat(devEdgLen,[nMembers+1 1 1]) ;
    hd.Seeds(crossSeedNumber).DataFields.MeanEdgRot = repmat(meanEdgRot,[nMembers+1 1 1]) ;
    hd.Seeds(crossSeedNumber).DataFields.DevEdgRot = repmat(devEdgRot,[nMembers+1 1 1]) ;
    
%%
    
% Macroscopic felds
    MacroL11 = hd.Seeds(macroSeedNumber).tri2nod*squeeze(hd.Seeds(macroSeedNumber).DataFields.L11) ;
    MacroL22 = hd.Seeds(macroSeedNumber).tri2nod*squeeze(hd.Seeds(macroSeedNumber).DataFields.L22) ;
    
% Display
    clf ;
    nodesIndices = 1:nNodes ; 6:9:nNodes ; [1:14 29:42] ; 4:7:nNodes ;
    xdata =  (1:hd.nFrames)' ...
            ... MacroL11(nodesIndices,:)' ...
            ... MacroL22(nodesIndices,:)' ...
            ;
    ydata =  reshape(permute((meanEdgRot(nodesIndices,:,:,:)*180/pi),[3 1 2 4]),hd.nFrames,[]) ... 
            ... reshape(permute((devEdgRot(nodesIndices,:,:,:)*180/pi),[3 1 2 4]),hd.nFrames,[]) ...
            ... (-MacroL11(nodesIndices,:)./MacroL22(nodesIndices,:))' ...
            ;
    pl = plot(xdata,ydata) ;
    lgd = legend(pl,arrayfun(@(i)num2str(i),nodesIndices,'uniformoutput',false)) ;
    lgd.EdgeColor = 'k' ;
    
    