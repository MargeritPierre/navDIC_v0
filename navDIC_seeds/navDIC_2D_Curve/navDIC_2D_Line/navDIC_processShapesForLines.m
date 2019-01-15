function H = navDIC_processShapesForDistMesh(H)

    % Is there geometries to process ?
        Points = [] ;
        if isempty(H.Geometries) || ~any([H.Geometries.isValid]) 
            pts = findobj(H.Axes,'tag','PointsPreview') ;
            delete(pts) ;
            return ; 
        end
        
    % Retrieve custom Data
        customData = H.uiContextMenuData ;
        % Convert to double
            for s = 1:length(customData) ;
                customData(s).N = str2num(customData(s).N) ;
            end
        % Set variable bar length
            N = [customData.N] ;
            
    % Recursive Distance function
        for s = 1:length(H.Shapes)
            if ~H.Geometries(s).isValid ; continue ; end
            pos = H.Geometries(s).Position ;
            Points(:,1,s) = interp1([0;1],pos(:,1),linspace(0,1,N(s))') ;
            Points(:,2,s) = interp1([0;1],pos(:,2),linspace(0,1,N(s))') ;
        end
        
    % Plot points
        pts = findobj(H.Axes,'tag','PointsPreview') ;
        if isempty(pts)
            pts = plot(H.Axes,NaN,NaN,'.b') ;
        end
        pts.XData = Points(2:end-1,1) ;
        pts.YData = Points(2:end-1,2) ;
        pts.MarkerSize = 20 ;
        pts.Tag = 'PointsPreview' ;
        pts.HitTest = 'off' ;
        uistack(findobj(H.Axes,'type','group'),'top') ;
        
    % Return the new Handle structure
        H.Points = Points ;
        
end
        