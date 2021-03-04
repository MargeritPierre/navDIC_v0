function [IDs,valid] = selectCameras(hd,SelectionMode)
    % Init
        IDs = [] ;
        valid = false ;
    % Are there any cameras ?
        if isempty(hd.Cameras) ; return ; end
    % Is there only one camera ?
        if length(hd.Cameras)==1 
            IDs  = 1 ; 
            valid = true ; 
            return ;
        end
    % Otherwise, choose a camera
    listCams = {hd.Cameras.Name} ;
    [IDs,valid] = listdlg('PromptString','Select Cameras :',...
                                'SelectionMode',SelectionMode,...
                                'initialValue',1,...
                                'ListString',listCams) ;
end