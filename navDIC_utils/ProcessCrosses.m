%% PROCESS CROSSES TO COMPUTE MEMBERS ANGLES, ELONGATIONS, ETC.

% Parameters
    crossSeedNumber = 2 ;
    macroSeedNumber = 1 ;
    nMembers = 4 ;
    
% Retrieve data
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
    devEdgLen = EdgLen-meanEdgLen ;
% Edge Rotations
    EdgRot = angle(Edg(:,:,:,:)./Edg(:,:,1,:)) ;
    meanEdgRot = mean(EdgRot,4) ;
    devEdgRot = mean(abs(EdgRot-meanEdgRot),4) ;
    
% Macroscopic felds
    MacroL11 = hd.Seeds(macroSeedNumber).tri2nod*squeeze(hd.Seeds(macroSeedNumber).DataFields.L11) ;
    MacroL22 = hd.Seeds(macroSeedNumber).tri2nod*squeeze(hd.Seeds(macroSeedNumber).DataFields.L22) ;
    
% Display
    clf ;
    nodesIndices = [1:14 29:42] ; 1:nNodes ; 4:7:nNodes ;
    xdata = ... (1:hd.nFrames)' ...
            ... MacroL11(nodesIndices,:)' ...
             MacroL22(nodesIndices,:)' ...
            ;
    ydata = ... reshape(permute((meanEdgRot(nodesIndices,:,:,:)*180/pi),[3 1 2 4]),hd.nFrames,[]) ... 
            ... reshape(permute((devEdgRot(nodesIndices,:,:,:)*180/pi),[3 1 2 4]),hd.nFrames,[]) ...
             (-MacroL11(nodesIndices,:)./MacroL22(nodesIndices,:))' ...
            ;
    pl = plot(xdata,ydata) ;
    lgd = legend(pl,arrayfun(@(i)num2str(i),nodesIndices,'uniformoutput',false)) ;
    
    