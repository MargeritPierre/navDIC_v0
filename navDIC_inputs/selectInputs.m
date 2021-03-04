function [IDs,valid] = selectInputs(hd,SelectionMode)
    % Init
        IDs = [] ;
        valid = false ;
    % Are there any cameras ?
        if isempty(hd.DAQInputs.Inputs) ; return ; end
    % Is there only one camera ?
        if length(hd.DAQInputs.Inputs)==1 
            IDs  = 1 ; 
            valid = true ; 
            return ;
        end
    % Otherwise, choose a camera
    listIns = {hd.DAQInputs.Inputs.Name} ;
    [IDs,valid] = listdlg('PromptString','Select Inputs :',...
                                'SelectionMode',SelectionMode,...
                                'initialValue',1,...
                                'ListString',listIns) ;
end