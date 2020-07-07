%% PROCESS CROSSES TO COMPUTE MEMBERS ANGLES, ELONGATIONS, ETC.

% Parameters
    crossSeedNumber = 3 ;
    macroSeedNumber = 2 ;
    unitCellsSeedNumber = 6 ;
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
    
% MACRO DATA
    
% Macroscopic fields
    macro = hd.Seeds(macroSeedNumber).computeDataFields(true) ; % Data Fields on nodes
% Unit Cell fields
    nUC = size(hd.Seeds(unitCellsSeedNumber).Elems,1) ;
    % Which cross is in which unit cell ?
        nodeInCell = zeros(nNodes,1) ;
        Xuc = pkg.data.dataAtIndices(hd.Seeds(unitCellsSeedNumber).Points,hd.Seeds(unitCellsSeedNumber).Elems) ;
        P = hd.Seeds(macroSeedNumber).Points ;
        for uc = 1:size(hd.Seeds(unitCellsSeedNumber).Elems,1)
            in = inpolygon(P(:,1),P(:,2),Xuc(uc,:,1),Xuc(uc,:,2)) ;
            nodeInCell(in) = uc ;
        end
    % Transfer matrix
        uc2cross = sparse(find(nodeInCell),nodeInCell(nodeInCell>0),true,nNodes,nUC) ;
    % Data Fields on elements
        unitcells = hd.Seeds(unitCellsSeedNumber).computeDataFields(false) ;
    % Data Fields on crosses
        for field = fieldnames(unitcells)'
            if isnumeric(unitcells.(field{1}))
                sz = size(unitcells.(field{1})) ;
                if sz(1)~=nUC ; continue ; end
                unitcells.(field{1}) = reshape(uc2cross*reshape(unitcells.(field{1}),nUC,[]),[nNodes sz(2:end)]) ;
                unitcells.(field{1})(nodeInCell<1,:) = NaN  ;
            end
        end
        
    
%% Display
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
    
    
%% GENERATE FIGURES FOR THE PAPER
    crossSeed = hd.Seeds(crossSeedNumber) ;
    nodesIndices = find(~ismember(nodeInCell,[0 2])) ; 1:nNodes ;
    fr = 1:hd.nFrames - 0 ;
    ucColors = hsv(nUC) ; jet(nUC) ; linspecer(nUC) ;

    E11 = unitcells.E11(nodesIndices,:).' ;
    E22 = unitcells.E22(nodesIndices,:).' ;
    E12 = unitcells.E12(nodesIndices,:).' ;
    rot = reshape(crossSeed.DataFields.dAtot(nodesIndices(:)+[0:nMembers-1]*nNodes,:),[],nMembers,hd.nFrames) ;
    meanRot = reshape(mean(rot,2),[],hd.nFrames).' ;
    devRot = reshape(sqrt(mean(abs(rot-mean(rot,2)).^2,2)),[],hd.nFrames).' ;
    %
    %meanRot = meanEdgRot(nodesIndices,:).' ; 
    %devRot = devEdgRot(nodesIndices,:).' ;
    
    clf
    ax = gca; ax.ColorOrder = ucColors(nodeInCell(nodesIndices),:) ;
    %pl = plot(meanRot(fr,:)) ; %xlabel '$image$' ; ylabel 'Mean Cross Rotation ($^\circ$)' ;
    pl = plot3(E12(fr,:)*100,180/pi*meanRot(fr,:),E22(fr,:)*100) ; xlabel '$E_{12} (\%)$' ; ylabel 'Mean Cross Rotation ($^\circ$)' ; zlabel '$E_{22} (\%)$' ;
    %pl = plot3(E22(fr,:)*100,180/pi*meanRot(fr,:),E12(fr,:)*100) ; xlabel '$E_{22} (\%)$' ; ylabel 'Mean Cross Rotation ($^\circ$)' ; zlabel '$E_{12} (\%)$' ;
    %pl = plot3(E12(fr,:)*100,180/pi*devRot(fr,:),E22(fr,:)*100) ; xlabel '$E_{12} (\%)$' ; ylabel 'Mean Member Deviation ($^\circ$)' ; zlabel '$E_{22} (\%)$' ;
    %pl = plot3(E22(fr,:)*100,180/pi*devRot(fr,:),E12(fr,:)*100) ; xlabel '$E_{22} (\%)$' ; ylabel 'Mean Member Deviation ($^\circ$)' ; zlabel '$E_{12} (\%)$' ;
    %pl = plot3(E22(fr,:)*100,E11(fr,:)*100,E12(fr,:)*100) ; xlabel '$E_{22} (\%)$' ; ylabel '$E_{11} (\%)$' ; zlabel '$E_{12} (\%)$' ;
    
    %set(ax,'xdir','reverse','ydir','reverse')
    
    labels = strtrim(num2cell(num2str((1:1000)'),2)) ;
    %lbl = strcat('C',labels(nodeInCell(nodesIndices))) ;
    lbl = strcat('N',labels(nodesIndices),'-C',labels(nodeInCell(nodesIndices))) ;
    [pl.DisplayName] = deal(lbl{:}) ;
    
    %%
    set(pl,'MarkerIndices',round(linspace(1,hd.nFrames,20))) ;
    set(pl,'MarkerSize',8) ;
    if 1
        set(gca,'sortmethod','childorder')
        set(pl,'color',[1 1 1]*0.6,'linewidth',1.5) ;
        colors = [1 0 0 ; 0 0 1 ; 0 1 0] ;
        markers = {'square','^','o','pentagram'} ;
        cells = [5 8] ;
        nodes = [41 42 50 51 ; 59 60 68 69] ;
        for cc = 1:numel(cells)
            cellNum = cells(cc) ;
            plCell = findobj(pl,'-regexp','DisplayName',['C' num2str(cellNum)]) ;
            set(plCell,'color',colors(cc,:))
            uistack(plCell,'top')
            for nn = 1:size(nodes,2)
                plNode = findobj(pl,'-regexp','DisplayName',['N' num2str(nodes(cc,nn))]) ;
                plNode.Marker = markers{nn} ;
            end
        end
    else
    end
    
%% 
    nodesIndices = find(nodeInCell) ; 1:nNodes ;

    E11 = unitcells.E11(nodesIndices,:).' ; % [nFrames nNodes]
    E22 = unitcells.E22(nodesIndices,:).' ; % [nFrames nNodes]
    E12 = unitcells.E12(nodesIndices,:).' ; % [nFrames nNodes]
    
    erot = EdgRot(nodesIndices,:,:,:) ; % [nNodes 1 nFrames nMembers]
    erot = permute(erot,[3 1 4 2]) ; % [nFrames nNodes nMembers]
    erot = erot(:,:,2) ; % [nFrames nNodes nMembers]
    
    clf
    ax = gca; %ax.ColorOrder = ucColors(nodeInCell(nodesIndices),:) ;
    pl = plot3(repmat(E12,[1 size(erot,3)])*100,180/pi*erot(:,:),repmat(E22,[1 size(erot,3)])*100) ; xlabel '$E_{12} (\%)$' ; ylabel 'Member Rotation ($^\circ$)' ; zlabel '$E_{22} (\%)$' ;
    %pl = plot3(E22*100,180/pi*meanRot,E12*100) ; xlabel '$E_{22} (\%)$' ; ylabel 'Mean Cross Rotation ($^\circ$)' ; zlabel '$E_{12} (\%)$' ;
    %pl = plot3(E12*100,180/pi*devRot,E22*100) ; xlabel '$E_{12} (\%)$' ; ylabel 'Mean Member Deviation ($^\circ$)' ; zlabel '$E_{22} (\%)$' ;
    %pl = plot3(E22*100,180/pi*devRot,E12*100) ; xlabel '$E_{22} (\%)$' ; ylabel 'Mean Member Deviation ($^\circ$)' ; zlabel '$E_{12} (\%)$' ;
    %pl = plot3(E11*100,E22*100,E12*100) ; xlabel '$E_{11} (\%)$' ; ylabel '$E_{22} (\%)$' ; zlabel '$E_{12} (\%)$' ;
   
    
    
    