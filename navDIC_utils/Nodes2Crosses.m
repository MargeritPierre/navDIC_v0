%% CONVERT A GRID OF POINTS TO CROSSES

%% CROSSES WITH 5 POINTS
%
%      2
%    / | \
%   3--5--1
%    \ | /
%      4
%

% Parameters
    global hd
    % New Seed
        fromSeedNumber = 5 ;
        newSeedNumber = 6 ;
        newSeedName = 'CrossesGrid' ;
    % Crosses
        L = 40 ; % Member length (in pixels)
        AngleMembers = [0 90 180 270]*pi/180 ; % Angle of the cross members
        AngleShift = repmat([1;-1;1;-1;1;-1;1;-1;1;-1;-1;1;-1;1;-1;1;-1;1;-1;1]*5*pi/180,1000,1) ; 60*pi/180 ;
    
% Create the new Seed
    hd.Seeds(newSeedNumber) = copy(hd.Seeds(fromSeedNumber)) ;
    hd.Seeds(newSeedNumber).Name = newSeedName ;

% Create crosses
    nNodes = size(hd.Seeds(newSeedNumber).Points,1) ;
    nMembers = numel(AngleMembers) ;
    % New nodes and displacements
        hd.Seeds(newSeedNumber).MovingPoints = repmat(hd.Seeds(newSeedNumber).MovingPoints,[nMembers+1 1 1]) ; 
        AngleShift = reshape(AngleShift(1:nNodes),[],1).*ones(nNodes,1) ;
    % Instantaneous Position shift
        for m = 1:nMembers % the last one will be the center node
            theta = AngleShift(:)+AngleMembers(m) ;
            % Translate 
                hd.Seeds(newSeedNumber).MovingPoints((1:nNodes)+(m-1)*nNodes,:,:) = L*[cos(theta) sin(theta)] + hd.Seeds(newSeedNumber).MovingPoints((1:nNodes)+(m-1)*nNodes,:,:) ;
        end
    % Other point data
        hd.Seeds(newSeedNumber).Points = hd.Seeds(newSeedNumber).MovingPoints(:,:,1) ; 
    % Elements
        hd.Seeds(newSeedNumber).Elems = [] ;
        for t = 1:nMembers
            if t==nMembers
                indP = [nMembers+1 t 1] ;
            else
                indP = [nMembers+1 t 1+t] ;
            end
            hd.Seeds(newSeedNumber).Elems((1:nNodes)+(t-1)*nNodes,:) = (1:nNodes)' + (indP(:)'-1)*nNodes ;
        end
        
% Other Data Fields
    hd.Seeds(newSeedNumber).computeDataFields ;
        
        
% DISPLAY THE RESULT
    cla ;
    imagesc(hd.Seeds(newSeedNumber).refImgs{1}) ;
    axis ij
    axis equal 
    axis tight
    axis off
    colormap gray
    pa = patch('Vertices',hd.Seeds(newSeedNumber).Points,'Faces',hd.Seeds(newSeedNumber).Triangles,'facecolor','w','facealpha',.1,'edgecolor','r') ;
    pts = plot(hd.Seeds(newSeedNumber).Points(:,1),hd.Seeds(newSeedNumber).Points(:,2),'.r','markersize',15) ;
        
        
        
        

%% CROSSES WITH 12 POINTS
%
%           3--2
%           |  |
%       5---4  1--12
%       |          |
%       6---7 10--11
%           |  |
%           8--9
%

% Parameters
    % New Seed
        fromSeedNumber = 1 ;
        newSeedNumber = 3 ;
        newSeedName = 'Crosses2' ;
    % Crosses
        L = 40 ; % Member length (in pixels)
        W = 10 ; % Member width (in pixels)
        AngleMembers = [0 90 180 270]*pi/180 ; % Angle of the cross members
        AngleShift = repmat([1;-1]*5*pi/180,1000,1) ; 60*pi/180 ;
    
%% Create the new Seed
    hd.Seeds(newSeedNumber) = copy(hd.Seeds(fromSeedNumber)) ;
    hd.Seeds(newSeedNumber).Name = newSeedName ;

% Create crosses
    nNodes = size(hd.Seeds(newSeedNumber).Points,1) ;
    nMembers = numel(AngleMembers) ;
    % New nodes and displacements
        hd.Seeds(newSeedNumber).MovingPoints = repmat(hd.Seeds(newSeedNumber).MovingPoints,[nMembers*3 1 1]) ; 
        AngleShift = reshape(AngleShift(1:nNodes),[],1).*ones(nNodes,1) ;
    % Instantaneous Position shift
        for m = 1:nMembers*4 % the last one will be the center node
            pos = []
            theta = AngleShift(:)+AngleMembers(m) ;
            % Translate 
                hd.Seeds(newSeedNumber).MovingPoints((1:nNodes)+(m-1)*nNodes,:,:) = L*[cos(theta) sin(theta)] + hd.Seeds(newSeedNumber).MovingPoints((1:nNodes)+(m-1)*nNodes,:,:) ;
        end
    % Other point data
        hd.Seeds(newSeedNumber).Points = hd.Seeds(newSeedNumber).MovingPoints(:,:,1) ; 
        hd.Seeds(newSeedNumber).Displacements = hd.Seeds(newSeedNumber).MovingPoints - hd.Seeds(newSeedNumber).Points ;
    % Triangles
        hd.Seeds(newSeedNumber).Triangles = [] ;
        for t = 1:nMembers
            if t==nMembers
                indP = [nMembers+1 t 1] ;
            else
                indP = [nMembers+1 t 1+t] ;
            end
            hd.Seeds(newSeedNumber).Triangles((1:nNodes)+(t-1)*nNodes,:) = (1:nNodes)' + (indP(:)'-1)*nNodes ;
        end
        
% Other Data Fields
    hd.Seeds(newSeedNumber).computeDataFields ;
        
        
% DISPLAY THE RESULT
    cla ;
    imagesc(hd.Seeds(newSeedNumber).refImgs{1}) ;
    axis ij
    axis equal 
    axis tight
    axis off
    colormap gray
    pa = patch('Vertices',hd.Seeds(newSeedNumber).Points,'Faces',hd.Seeds(newSeedNumber).Triangles,'facecolor','w','facealpha',.1,'edgecolor','r') ;
    pts = plot(hd.Seeds(newSeedNumber).Points(:,1),hd.Seeds(newSeedNumber).Points(:,2),'.r','markersize',15) ;
        
        
        
        