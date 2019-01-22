function [IDs,valid] = selectSeeds(hd,SelectionMode)
    % Init
        IDs = [] ;
        valid = false ;
    % Are there any seeds ?
        if isempty(hd.Seeds) ; return ; end
    % Is there only one seed ?
        if length(hd.Seeds)==1 
            IDs  = 1 ; 
            valid = true ; 
            return ;
        end
    % Otherwise, choose a seed
        listSeeds = {hd.Seeds.Name} ;
        [IDs,valid] = listdlg('PromptString','Select a Seed :',...
                                    'SelectionMode',SelectionMode,...
                                    'initialValue',1,...
                                    'ListString',listSeeds) ;
end