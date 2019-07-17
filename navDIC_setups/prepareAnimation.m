function out = prepareAnimation(hd)

    disp '------ ANIMATION -------'
    
    % Initialize
        out.Valid = false ;
    
    % Select the preview to be animated
        [IDs,valid] = selectPreviews(hd,'multiple') ;
        if ~valid ; warning('NO VALID PREVIEW TO EXPORT') ; return ; end
        
    % Retrieve the preview
        prev = [hd.Previews{IDs}] ;
        if any(~[prev.isValid]) ; warning('AT LEAST ONE OF THE SELECTED PREVIEWS IS NOT VALID') ; return ; end
        
    % Set the previews as current figures
        figs = [prev.fig] ;
        for fi = figs
            figure(fi) ;
        end
        
    % Ask for the file to be saved
        defaultName = [hd.WorkDir.Path,'\','Video','_',char(regexprep(string(datetime),{' ','-',':'},'_')),'.avi'] ;
        [file,path] = uiputfile('*.avi','WHERE TO SAVE THE ANIMATION ?',defaultName) ;
        if file==0 ; warning('THE SELECTED FILE IS NOT VALID') ; return ; end
        VideoFile = [path file] ;
        
    % Try to choose the right image stacking method (multi-previews)
        FR = {} ;
        sizes = [] ;
        for p = 1:length(figs)
            FR{end+1} = getframe(figs(p)) ;
            FR{end} = FR{end}.cdata ;
            sizes(end+1,:) = size(FR{end}(:,:,1)) ;
        end
        horizontalStackRatio = max(sizes(:,1))/sum(sizes(:,2)) ;
        verticalStackRatio = sum(sizes(:,1))/max(sizes(:,2)) ;
        if abs(horizontalStackRatio-1)<abs(verticalStackRatio-1)
            stack = 'horizontal' ;
        else
            stack = 'vertical' ;
        end
        
    % Ask for other parameters
        % Parameters
            prompt = {'Quality (%):','Frames Recorded:'};
            definput = {'75',['[1:1:',num2str(hd.nFrames),']']};
            if numel(figs)>1
                prompt{end+1} = 'Stacking (''horizontal'' or ''vertical'')' ;
                definput{end+1} = stack ;
            end
        % Ask
            dims = [1 35];
            answer = inputdlg(prompt,'Animation Parameters',dims,definput) ;
            if isempty(answer) ; return ; end
        % Fill Values
            params.VideoQuality = str2num(answer{1}) ;
            FramesRecorded = str2num(answer{2}) ;
            if numel(figs)>1
                stack = lower(answer{3});
            end
        
    % Other parameters (has to be putted in an interface later...)
        
    
    % Create the Video Writer
        writerObj = VideoWriter(VideoFile) ;
        writerObj.Quality = params.VideoQuality ;
        
    % Frame taking function
        figPositions = cat(1,figs.Position) ;
        switch stack
            case 'horizontal'
                [~,stackOrder] = sort(figPositions(:,1),'ascend') ;
            case 'vertical'
                [~,stackOrder] = sort(figPositions(:,2),'descend') ;
        end
        IDs = IDs(stackOrder) ;
        figs = figs(stackOrder) ;
        sizes = sizes(stackOrder,:) ;
        
    % Output structure
        out.Valid = true ;
        out.prevID = IDs ;
        out.figs = figs ;
        out.writerObj = writerObj ;
        out.params = params ;
        out.FramesRecorded = FramesRecorded ;
        out.sizes = sizes ;
        out.stack = stack ;
            
end