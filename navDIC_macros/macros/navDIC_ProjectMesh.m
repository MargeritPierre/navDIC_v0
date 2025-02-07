classdef navDIC_ProjectMesh < navDIC_AbstractMacro
%NAVDIC_PROJECTMESH project a master seed to a slave seed
properties
    MasterSeed navDICSeed = navDICSeed.empty
    SlaveSeed navDICSeed = navDICSeed.empty
    P % the projection matrix
end

methods
    function this = navDIC_ProjectMesh()
    % Class constructor
        if isempty(this.Name) ; this.Name = regexprep(class(this),'navDIC_','') ; end
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        seedNames = strjoin(strcat(string(num2cell(1:numel(hd.Seeds))),':',{hd.Seeds.Name}),newline) ;
        [~,MasterID] = ismember(this.MasterSeed,hd.Seeds) ;
        [~,SlaveID] = ismember(this.SlaveSeed,hd.Seeds) ;
        seedIDs = inputdlg(...
                    {"Master Seed:" + newline + seedNames ; "Slave Seed:" + newline + seedNames ; } ...
                    ,'Seed projection parameters' ...
                    ,1 ...
                    ,{num2str(MasterID) ; num2str(SlaveID)} ...
                    ) ;
        MasterID = str2double(seedIDs{1}) ;
        SlaveID = str2double(seedIDs{2}) ;
        if any(~ismember([MasterID SlaveID],1:numel(hd.Seeds))) ; return ; end
        this.MasterSeed = hd.Seeds(MasterID) ;
        this.SlaveSeed = hd.Seeds(SlaveID) ;
        this.P = this.MasterSeed.interpMat(this.SlaveSeed.Points) ;
    end

    function hd = run(this,hd)
    % Run the macro
        fr = hd.CurrentFrame ;
        um = this.MasterSeed.MovingPoints(:,:,fr)-this.MasterSeed.Points ;
        xs = this.SlaveSeed.Points ;
        this.SlaveSeed.MovingPoints(:,:,fr) = xs + this.P*um ;
        this.SlaveSeed.computeDataFields([],fr) ;
    end
end

end

