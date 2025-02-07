classdef (Abstract) navDIC_DICMacro < navDIC_AbstractMacro
%NAVDIC_DICMacro Base class for navDIC DIC macros

properties
    % User-settable properties
    Seed navDICSeed = navDICSeed.empty % the corresponding seed
    RefFrame = 0 % reference frame number. if==0, tak the seed's ref frame
    WeightPreviousImage double = 0 % refImg = (1-w)*refImg + w*prevImg
    GaussianFilterSize double = 0 ; % image gaussian pre-filter
    InitWithExistingData logical = true % initialize with the current seed's .MovingPoints
    UsePreviousVelocity logical = true % add the previous motion to the initial guess
    DispComp = 'both' ; % Estimated displacement components
    % Other properties
    RefImgs cell % the reference image(s)
end

methods
    function this = navDIC_DICMacro()
    % Class constructor
        if isempty(this.Name) ; this.Name = regexprep(class(this),'navDIC_','') ; end
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
    % First, choose the seed
        if isempty(hd.Seeds) 
            errordlg('A DIC Seed needs to be created first !','ERROR') ; 
            return ; 
        end
        seedID = inputdlg(...
                    strjoin(strcat(string(num2cell(1:numel(hd.Seeds))),':',{hd.Seeds.Name}),newline) ...
                    ,'Choose a seed on which to perform DIC' ...
                    ,1 ...
                    ,{num2str(find(ismember(hd.Seeds,this.Seed)))} ...
                    ) ;
        seedID = str2double([seedID{:}]) ;
        if ~any(ismember(seedID,1:numel(hd.Seeds))) ; return ; end
    % Retrieve the seed parameters if changed
        if all(this.Seed~=hd.Seeds(seedID))
            this.Seed = hd.Seeds(seedID) ;
            this.Name = [regexprep(class(this),'navDIC_','') '_' this.Seed.Name] ;
            this.Enable = this.Seed.compDisp ;
            switch this.Seed.displMode
                case 'rel' ; this.WeightPreviousImage = 1 ;
                otherwise ; this.WeightPreviousImage = 0 ;
            end
            this.InitWithExistingData = this.Seed.useExistingDisp ;
            this.UsePreviousVelocity = this.Seed.useExistingDisp ;
        end
    % Then set the DIC parameters
        defInputs = { ...
                        'Name' , this.Name ...
                        ; 'Enable [0/1]' , num2str(this.Enable) ...
                        ; 'Reference frame (0 to take the seed''s reference)' , num2str(this.RefFrame) ...
                        ; ['Weight of the previous image in the reference:' newline 'refImg = (1-w)*refImg + w*prevImg'] , num2str(this.WeightPreviousImage) ...
                        ; 'Gaussian filter size' , num2str(this.GaussianFilterSize) ...
                        ; 'Initialize with existing data [0/1]' , num2str(this.InitWithExistingData) ...
                        ; 'Use the previous velocity [0/1]' , num2str(this.UsePreviousVelocity) ...
                        ; 'Estimated displacement components [both/X/Y]' , this.DispComp ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'Common DIC parameters' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
        if isempty(out) ; return ; end
        this.Name = out{1} ;
        this.Enable = str2num(out{2}) ;
        switch str2num(out{3})
            case 0 ; this.RefImgs = this.Seed.refImgs ;
            otherwise ; this.RefImgs = hd.Images{this.Seed.CamIDs}(str2num(out{3})) ;
        end
        this.WeightPreviousImage = str2num(out{4}) ;
        this.GaussianFilterSize = str2num(out{5}) ;
        this.InitWithExistingData = str2num(out{6}) ;
        this.UsePreviousVelocity = str2num(out{7}) ;
        this.DispComp = out{8} ;
    % Filter the reference image
        this.RefImgs = this.processImgs(this.RefImgs);
    end
end

methods (Abstract)
% THESE METHODS NEED TO BE INSTANCIED BY SUBCLASSES
    X = updateDIC(this,X,img,IMG) % Update the configuration X using DIC performed from current images "img" to references "IMG"
    I = transformImage(this,I,x,X) % Transform the image(s) I from configuration X to x
end

methods (Sealed)
    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        if ~this.Enable ; return ; end
    % Reference and current frames
        if this.RefFrame==0 ; RF = this.Seed.RefFrame ; 
        else ; RF = this.RefFrame ; end
        CF = hd.CurrentFrame ;
        if CF<this.Seed.RefFrame ; return ; end
    % Add dummy displacement data if needed
        if size(this.Seed.MovingPoints,3)<CF
            this.Seed.MovingPoints(:,:,hd.CurrentFrame) = NaN ;
        end
    % Initialize the configuration
        if CF==RF % with the reference config
            this.Seed.MovingPoints(:,:,CF) = this.Seed.Points ;
            this.Seed.computeDataFields ; % initalize datafields
            return ; % the reference config is known
        else % with the previous config
            X0 = this.Seed.MovingPoints(:,:,CF-1) ;
        end
    % Take the available data if needed
        dataAvailable = any(~isnan(this.Seed.MovingPoints(:,:,CF)),[1 2]) ;
        useData = this.InitWithExistingData && dataAvailable ;
        if useData
            X0 = this.Seed.MovingPoints(:,:,CF) ;
        end
    % Add previous velocity
        if this.UsePreviousVelocity ... if asked for
            && CF-RF>=2 ...if there is at least two previous known configs
            && ~useData % if there is no other data that has been used
            X0 = X0 + diff(this.Seed.MovingPoints(:,:,CF+(-2:-1)),1,3) ;
        end
    % Get and process the current image(s)
        imgs = hd.Images{this.Seed.CamIDs}(CF) ;
        imgs = this.processImgs(imgs) ;
    % Update the DIC
        X = this.updateDIC(X0,imgs,this.RefImgs) ;
    % Update the reference images if needed
        if this.WeightPreviousImage
            imgs = this.transformImage(imgs,X,X0) ;
            for ii = 1:numel(this.RefImgs)
                this.RefImgs{ii} = this.RefImgs{ii}*(1-this.WeightPreviousImage) + imgs{ii}*this.WeightPreviousImage ;
            end
        end
    % Set the seed's moving points
        this.Seed.MovingPoints(:,:,CF) = X ;
        this.Seed.computeDataFields([],CF) ;
    end
end

methods
    function imgs = processImgs(this,imgs)
    % Process images before DIC
        isCell = iscell(imgs) ;
    % Convert to double and normalize
        for ii = 1:numel(imgs)
            imgs{ii} = double(imgs{ii})*(1/max(getrangefromclass(imgs{ii}(1)))) ;
        end
    % Convert to black and white
        for ii = 1:numel(imgs)
            imgs{ii} = mean(imgs{ii},3) ;
        end
    % Apply a gaussian filter on images
        if all(this.GaussianFilterSize<=1) ; return ; end
        sz = this.GaussianFilterSize.*[1 1] ;
    % ND Gaussian kernel
        sig = (sz+1)./log(5*sz) ; % optimized variance
        xx = arrayfun(@colon,-sz,sz,'UniformOutput',false) ;
        [xx{:}] = ndgrid(xx{:}) ;
        xx = cat(numel(sz)+1,xx{:})./reshape(sig,[ones(1,numel(sz)) numel(sz)]) ;
        kern = exp(-sum(xx.^2,numel(sz)+1)) ; 
    % Convolve
        for ii = 1:numel(imgs)
            imgs{ii} = convn(imgs{ii},kern,'same') ;
        end
    end
end

%% IMAGE INTERPOLATION
methods
    function val = bilinearInterp(this,img,xx,extrapVal)
    % Bilinear interpolation of an image at real-valued coordinates xx
    % faster than gg = interp2(img,xx(:,1),xx(:,2),'linear') ;
    % xx = [nValues 2] = [jj(:) ii(:)] ;
    % with ii and jj resp. the (real-valued) row and column indices
    % extrapVal: extrapolation value (default: NaN)

        if nargin<4 ; extrapVal = NaN ; end

        % Image informations
        [nI,nJ,~] = size(img) ;

        % Valid coordinates
        valid = xx(:,1)<=nJ ...
                & xx(:,1)>=1 ...
                & xx(:,2)<=nI ...
                & xx(:,2)>=1 ;

        % Dummy values
        xx(~valid,:) = 1 ;

        % Integer part of the cordinates
        ji = floor(xx) ;
        ji(ji(:,1)==nJ,1) = nJ-1 ;
        ji(ji(:,2)==nI,2) = nI-1 ;

        % Neightboring pixels
        p1 = ji(:,2)+nI*(ji(:,1)-1) ;
        p2 = p1 + nI ; 
        p3 = p2+1 ; 
        p4 = p1+1 ;

        % residual coordinates
        dx = xx-ji ;

        % bilinear interpolation
        val = img(p1).*(1-dx(:,1)).*(1-dx(:,2)) ...
            + img(p2).*dx(:,1).*(1-dx(:,2)) ...
            + img(p3).*dx(:,1).*dx(:,2) ...
            + img(p4).*(1-dx(:,1)).*dx(:,2) ;

        % Non-valid values
        val(~valid) = extrapVal ;
    end
end


end

