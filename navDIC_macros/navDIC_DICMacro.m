classdef (Abstract) navDIC_DICMacro < navDIC_AbstractMacro
%NAVDIC_DICMacro Base class for navDIC DIC macros

properties
    % User-settable properties
    Seed navDICSeed = navDICSeed.empty
    RefFrame = 0
    WeightPreviousImage double = 0
    InitWithExistingData logical = true
    UsePreviousVelocity logical = true
    GaussianFilterSize double = 0 ;
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
    % Filter the reference image
        this.RefImgs = this.processImgs(this.RefImgs);
    end
end

methods (Abstract)
% THESE METHODS NEEDS TO BE INSTANCIED BY SUBCLASSES
    X = updateDIC(this,X,imgs) % Update the configuration X using DIC performed on imgs
    I = transformImage(this,X,x,I) % Transform an image(s) I from configuration X to x
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
        X = updateDIC(this,X0,imgs) ;
    % Update the reference images if needed
        if this.WeightPreviousImage
            imgs = transformImage(this,X0,X,imgs) ;
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

end

