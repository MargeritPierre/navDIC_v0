function [IDs,valid,figIDs] = selectPreviews(hd,SelectionMode)
    % Init
        IDs = [] ;
        valid = false ;
    % Are there any previews ?
        if isempty(hd.Previews) ; return ; end
    % For previews with multiple figures
        nFig = cellfun(@numel,{hd.Previews.fig}) ; % number of figures in each preview
        prevID = repelem(1:numel(hd.Previews),1,nFig); % preview ID
        figID = arrayfun(@colon,1+nFig*0,nFig,'UniformOutput',false) ;
        figID = cat(2,figID{:}) ; % figure ID
    % Is there only one preview figure ?
        if sum(nFig)==1 
            IDs  = 1 ; 
            valid = true ; 
            return ;
        end
    % Otherwise, choose preview figure(s)
        listPreviews = [hd.Previews.fig] ;
        listPreviews = {listPreviews.Name} ;
        [IDs,valid] = listdlg('PromptString','Select a Preview :',...
                                    'SelectionMode',SelectionMode,...
                                    'initialValue',1,...
                                    'ListString',listPreviews) ;
    % Return previews and figures
        IDs = prevID(IDs) ;
        figIDs = figID(IDs) ;
end