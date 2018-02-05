%% IMPLEMENTS DIGITAL IMAGE CORRELATION ON PREVIOUSLY ACQUIRED IMAGES
clc
close all
clear all

%% LOAD DATA FOR DIC
    clc
    close all
    clear all

% Working Directory
    [path] = uigetdir('SELECT THE DIRECTORY OF IMAGES') ;
    if path==0 ; return ; end
    
% Image files present in the directory (in the format: commonName_id.ext)
    % Possible image extensions
        imgExt = {'tif','png','bmp','jpg','jpeg'} ;
    % Get files
        files = dir('*.anImpossibleExtension') ;
        for i = 1:length(imgExt)
            f = dir([path,'/*.',imgExt{i}]) ;
            files(end+(1:length(f))) = f ;
        end
    % Keep file names only
        fileNames = {files.name} ;
    % Get the common name and extension
        str = strsplit(fileNames{1},'_') ;
        if length(str)==1
            commonName = str{1} ;
        else
            commonName = strjoin(str(1:end-1),'_') ;
        end
        [~,~,ext] = fileparts(str{end}) ;
    % Get image ids
        idSTR = {} ;
        idNUM = [] ;
        for i = 1:length(fileNames)
            idSTR{i} = fileNames{i}(length(commonName)+2:end-length(ext)) ;
            if ~isempty(str2num(idSTR{i}))
                idNUM(i) = str2num(idSTR{i}) ;
            else
                idNUM(i) = NaN ;
            end
        end
        [idNUM,ind] = sort(idNUM(~isnan(idNUM))) ;
        idSTR = idSTR(ind(~isnan(idNUM))) ;
        nImgs = length(idSTR) ;
        
% Load images
    imgs = {} ;
    for i = 1:nImgs
        imgs{i} = imread([path,'/',commonName,'_',idSTR{i},ext]) ;
    end
    
% Load ROI
    ROI = imread([path,'/',commonName,'_ROI',ext]) ;
    ROI = uint8(logical(ROI(:,:,1))) ;
    
% Load TimeLine
    try
        load([path,'/',commonName,'_time.mat']) ;
    end
    
%% (OPTIONAL !!!!!) SET REGION OF INTEREST
    % Take a picture as refImg for ROI setting
        roiImg = imgs{1};
    % Open the GUI for ROI setting
        ROI = SetROI(roiImg,0) ;
    % Save the ROI as image
%         fileNameROI = [path,file,'_ROI',ext] ;
%         imwrite(ROI,fileNameROI);
    

%% PERFORM DIC
    
% USER PARAMETERS
    d = 5 ;
    CorrSize = 20 ; % Imagette size in pixels
    disp_Mode = 'abs' ; % 'abs' or 'rel'
    width_pln_strains = 8*d ; % max distance of ctrl pts to compute strains
    % Control Points parameters
        dx = d ;
        dy = d ;

% First grid
    im0 = imgs{1} ;
    [ly,lx] = size(im0) ;
    [X0,Y0] = meshgrid(1:dx:lx,1:dy:ly) ;
    % center the grid
        X0 = round(X0-mean(X0(:))+lx/2) ;
        Y0 = round(Y0-mean(Y0(:))+ly/2) ;

% Valid points inside ROI
    valid = logical(ROI(Y0(:,1)',X0(1,:))) ;
    X0 = X0(valid) ;
    Y0 = Y0(valid) ;

% Used grid
    Pts0 = [X0 Y0] ;
    nPts = size(Pts0,1) ;
    display([num2str(nPts),' correlation points'])
    
% Build delaunay mesh
    tri = delaunay(X0,Y0) ;
    % Delete triangles with long edges
        edg = sqrt((X0(tri(:,[2:end,1]))-X0(tri)).^2+(Y0(tri(:,[2:end,1]))-Y0(tri)).^2) ;
        OK = sum(edg<=norm([dx,dy])*1.01,2)==3 ;
        tri = tri(OK,:) ;

% Figure
    m = 0.02 ;
    fig = figure('windowstyle','docked') ;
        axImg = axes('position',[m m 1-2*m 1-2*m]) ;
            axImg.YDir = 'reverse' ;
            axis equal
            axis tight
            axis off
            im = imagesc(im0) ;
            im.CData = repmat(im0,[1 1 3]) ;
            pl = plot(X0,Y0,'.r','markersize',8) ;
            srf = trisurf(tri,X0,Y0,0*Y0,'facealpha',0.4,'facecolor','interp','edgecolor','none') ;
            clrbr = colorbar('location','manual','position',axImg.Position.*[0 1 0 1]+[.82 0 0.05 0]) ;
        drawnow ;
        pause(.01) ;

%% Perform DIC
    delete(pl) ; clear('pl') ;
    Pts = Pts0 ;
    U = zeros([size(Pts),nImgs]) ;
    E = zeros([nPts 3 nImgs]) ; % Exx,Eyy,Exy
    for i = 2:1:nImgs
        display(['Img ',num2str(i),'/',num2str(nImgs)])
        % Get Valid Points
            valid = ~any(isnan(Pts),2) ;
        % Perform DIC
            t = tic ;
            switch disp_Mode
                case 'abs'
                    Pts(valid,:) = my_cpcorr(round(Pts(valid,:)),Pts0(valid,:),imgs{i},imgs{1},CorrSize) ;
                case 'rel'
                    Pts(valid,:) = my_cpcorr(Pts(valid,:),Pts(valid,:),imgs{i},imgs{i-1},CorrSize) ;
            end
            toc(t) 
        % Compute Displacements
            U(:,:,i) = (Pts-Pts0) ;
        % Compute Strains
            for p = 1:nPts
                % Init
                    E(p,:,i) = NaN ;
                    if isnan(U(p,1,i)) ; continue ; end ;
                % Get neightoring Pts
                    pt = Pts(p,:) ;
                    r = sqrt(sum((Pts-repmat(pt,[nPts 1])).^2,2)) ;
                    indFit = r<width_pln_strains ;
                    nfit = size(indFit(indFit),1) ;
                    if nfit<3 ; continue ; end ;
                % Fit a plane
                    %fitPts = Pts(indFit,:) ;
                    xy = Pts(indFit,:)-repmat(pt,[nfit 1]) ;
                    % Plane : a+bx+cy = Ux ; d+ex+fy = Uy ; 
                        A = [xy(:,1) xy(:,2) ones(nfit,1)] ;
                        bx = U(indFit,1,i) ;
                        by = U(indFit,2,i) ;
                        Px = A\bx ;
                        Py = A\by ;
                % Strains
                    E(p,1,i) = Px(1) ;
                    E(p,2,i) = Py(2) ;
                    E(p,3,i) = .5*(Py(1)+Px(2)) ;
            end
            toc(t)
        % Update Plots
            if exist('pl')
                pl.XData = Pts(:,1) ;
                pl.YData = Pts(:,2) ;
            end
            im.CData = repmat(imgs{i},[1 1 3]) ;
            srf.Vertices = [Pts zeros(nPts,1)] ;
            srf.FaceVertexCData = E(:,2,i)*100 ; sqrt(sum(U.^2,2)) ;
            drawnow ;
            pause(.01) ;
    end
    
%% PLOT OTHER DATA
    srf.FaceVertexCData = E(:,1,end)*100
    
    
%% VALEURS MOYENNES
    % Exx
        Exx = E(:,1,:) ;
        Exx(isnan(Exx)) = 1i ;
        Exx = sum(real(Exx),1)./sum(imag(Exx)==0,1) ;
        Exx = squeeze(Exx) ; 
    % Exy
        Exy = E(:,3,:) ;
        Exy(isnan(Exy)) = 1i ;
        Exy = sum(real(Exy),1)./sum(imag(Exy)==0,1) ;
        Exy = squeeze(Exy) ;
    % Eyy
        Eyy = E(:,2,:) ;
        Eyy(isnan(Eyy)) = 1i ;
        Eyy = sum(real(Eyy),1)./sum(imag(Eyy)==0,1) ;
        Eyy = squeeze(Eyy) ;
    % FORCE
        FORCE = ones(size(Exx)) ;
        try
            FORCE = load([path,'/',commonName,'_inputs.mat']) ;
            FORCE = FORCE.data ;
        end
    
    figure ;
    plot(time,Exx) ;
    plot(time,Eyy) ;
    plot(time,Exy) ;
    legend('Exx','Eyy','$\epsilon_{xy}$')
    
    
%% POISSON AND YOUNG
    i = [1:nImgs]
    figure ;
    plot(time(i),-Exx(i)./Eyy(i))
    
    figure ;
    S = 11.4*0.95 ;
    plot(Eyy(i),FORCE(i)/S)
        














