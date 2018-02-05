function [ROI,Geometries] = SetROI(refImg,varargin)
% varargin is used in cases where no base rectangle is needed.

    global handles 
    handles = {} ;
    handles.Shapes = {} ;
    handles.Geometries = [] ;
    handles.RefImg = refImg ;
    handles.ROI = refImg(:,:,1).*0+1 ;
    
    % Basis rectangle ?
        if isempty(varargin) 
            handles.BasisRect = true ;
        else
            handles.BasisRect = false ;
        end
    
    initFig() ;
    updateROI() ;
    
    while ishandle(handles.Figure) 
        drawnow ;
        pause(.05) ;
    end
    
    ROI = handles.ROI ;
    Geometries = handles.Geometries ;
end


function initFig
    global handles 
    
    marg = 5 ; % Pixels
    refImg = handles.RefImg ;
    
    % Figure
        handles.Figure = figure('tag','figROI','Name','SETTING THE R.O.I') ;%,'toolbar','none','windowstyle','docked') ;
            maxSize = get(0,'defaultfigureposition') ; maxSize = maxSize(3:4) ;
            Ratio = max((size(refImg')+2*marg)./(maxSize-[0 55])) ;
            handles.Figure.Position([4,3]) = round(size(refImg)/Ratio)+2*marg ;
        
    % Image    
        handles.Axes = axes('position',[0 0 1 1]) ;
            handles.Axes.Units = 'pixels' ;
            handles.Axes.Position = handles.Axes.Position + [marg marg -2*marg -2*marg] ;
            handles.Axes.Units = 'normalized' ;
            handles.Axes.YDir = 'reverse' ;
            handles.Img = imagesc(refImg) ;
            handles.Img.CData = repmat(refImg,[1 1 3]) ;
            axis equal
            axis tight
            axis off
            
    % ROI
        handles.SrfROI = imagesc(refImg) ;
        handles.SrfROI.CData = repmat(handles.ROI,[1 1 3]).*0+255 ;
            
    % Basis Rectangle if needed
        if handles.BasisRect
            newShape('Rectangle','+',[0 0 size(handles.Img.CData,2) size(handles.Img.CData,1)]) ;
            handles.Shapes{end}.Deletable = 0 ;
        end
            
    % ToolBar    
        handles.ToolBar = uitoolbar(handles.Figure) ;
        addButtons ;  
end


function newShape(shape,msg,position)
    global handles
    
    switch shape
        case 'Rectangle'
            if isempty(position)
                handles.Shapes{end+1} = imrect(handles.Axes) ;
            else
                handles.Shapes{end+1} = imrect(handles.Axes,position) ; 
            end
            fcn = makeConstrainToRectFcn('imrect',handles.Axes.XLim, handles.Axes.YLim) ;
        case 'Ellipse'
            handles.Shapes{end+1} = imellipse(handles.Axes) ;
            fcn = makeConstrainToRectFcn('imellipse',handles.Axes.XLim, handles.Axes.YLim) ;
        case 'Polygon'
            handles.Shapes{end+1} = impoly(handles.Axes) ;
            fcn = makeConstrainToRectFcn('impoly',handles.Axes.XLim, handles.Axes.YLim) ;
    end
    
    switch msg
        case '+'
            setColor(handles.Shapes{end},'b') ;
        case '-'
            setColor(handles.Shapes{end},'r') ;
    end
    set(handles.Shapes{end},'DisplayName',msg) ;
    setPositionConstraintFcn(handles.Shapes{end},fcn) ;
    addNewPositionCallback(handles.Shapes{end},@(pos)updateROI) ;
    updateROI ;
end


function updateROI
    global handles
    
    ROI = logical(handles.RefImg(:,:,1).*0) ;
    shapes = handles.Shapes ;
    for s = 1:length(shapes)
        try
            msg = get(shapes{s},'DisplayName') ;
            switch msg
                case '+'
                    ROI = ROI | logical(createMask(shapes{s},handles.SrfROI)) ;
                case '-'
                    ROI = ROI & logical((1-createMask(shapes{s},handles.SrfROI))) ;
                otherwise
            end
        end
    end
    handles.SrfROI.AlphaData = (1-ROI)*.5 ;
    handles.ROI = double(ROI) ;
    
    % Backup Geometries
        handles.Geometries = [] ;
        for s = 1:length(handles.Shapes)
            if ~isvalid(handles.Shapes{s}) ; continue ; end ;
            handles.Geometries(end+1).Class = class(handles.Shapes{s}) ;
            handles.Geometries(end).Bool = get(handles.Shapes{s},'DisplayName') ;
            handles.Geometries(end).Position = getPosition(handles.Shapes{s}) ;
        end
end


function addButtons
    global handles
    
    % Add Rectangle    
        handles.AddRect = uipushtool(handles.ToolBar) ;
        handles.AddRect.TooltipString = 'Add a Rectangle' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconAddRect.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.AddRect.CData = icon ;
        % Declare Callback
            handles.AddRect.ClickedCallback = @(src,evt)newShape('Rectangle','+',[]) ;
    % Add Rectangle    
        handles.RemRect = uipushtool(handles.ToolBar) ;
        handles.RemRect.TooltipString = 'Remove a Rectangle' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconRemRect.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.RemRect.CData = icon ;
        % Declare Callback
            handles.RemRect.ClickedCallback = @(src,evt)newShape('Rectangle','-',[]) ;
    
    % Add Ellipse    
        handles.AddElli = uipushtool(handles.ToolBar) ;
        handles.AddElli.TooltipString = 'Add an Ellipse' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconAddElli.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.AddElli.CData = icon ;
        % Declare Callback
            handles.AddElli.ClickedCallback = @(src,evt)newShape('Ellipse','+',[]) ;
    % Remove Ellipse    
        handles.RemElli = uipushtool(handles.ToolBar) ;
        handles.RemElli.TooltipString = 'Remove an Ellipse' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconRemElli.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.RemElli.CData = icon ;
        % Declare Callback
            handles.RemElli.ClickedCallback = @(src,evt)newShape('Ellipse','-',[]) ;
    
    % Add Polygon   
        handles.AddPoly = uipushtool(handles.ToolBar) ;
        handles.AddPoly.TooltipString = 'Add a Polygon' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconAddPoly.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.AddPoly.CData = icon ;
        % Declare Callback
            handles.AddPoly.ClickedCallback = @(src,evt)newShape('Polygon','+',[]) ;
    % Remove Polygon   
        handles.RemPoly = uipushtool(handles.ToolBar) ;
        handles.RemPoly.TooltipString = 'Remove a Polygon' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconRemPoly.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.RemPoly.CData = icon ;
        % Declare Callback
            handles.RemPoly.ClickedCallback = @(src,evt)newShape('Polygon','-',[]) ;
            
            
    % Validate ROI 
        handles.roiOK = uipushtool(handles.ToolBar) ;
        handles.roiOK.TooltipString = 'Save ROI and Close' ;
        % Draw Icon
            [icon,~,alpha] = imread('iconRoiOK.png') ; 
            icon = double(icon)/255 ;
            icon(repmat(alpha,[1 1 3])==0) = NaN ;
            handles.roiOK.CData = icon ;
        % Declare Callback
            handles.roiOK.ClickedCallback = @(src,evt)close(handles.Figure) ;
          
end
