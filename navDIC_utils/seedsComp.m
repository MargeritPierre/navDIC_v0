function Comp = seedsComp(seedHP, hd ) 

CamIDs = seedHP.CamIDs ;
nbImgs = size(seedHP.Strains,3);
% Raideur fil
raideurFil = input( 'Entree la raideur du fil (N) : ' ) ;

% longueur reel mortier
rep = inputdlg('utiliser image de calibration differente (oui/non)') ;
if strcmpi(rep{1},'oui') 
    [file,path] = uigetfile('*.tif',['selectionner l''image de calibration de la camera : ',...
        hd.Cameras(num).Name,' ? utiliser. ']) ;
    img = imread([path,file]) ;
else
    if iscell(hd.Images{1}{CamIDs})
        img = hd.Images{1}{CamIDs}{1} ;
    else
        img = hd.Images{1}{CamIDs} ;
    end
end


drawToolH = drawingTool('drawROI',true ...
       ,'background', img,'title','DrawingTool : Define the size of the mortar prism') ;   
if strcmpi(drawToolH.Geometries(1).Class, 'impoint')
   for p = 1:length(drawToolH.Geometries)
        pts(p,:) = obj.drawToolH.Geometries(p).Position ;
   end
elseif strcmpi(drawToolH.Geometries(1).Class, 'imline')
   pts = drawToolH.Geometries(1).Position(:,:) ;
end
Comp.CamIDs = CamIDs ;
Comp.seedName = seedHP.Name ;
Comp.Lmortar = sqrt(sum(diff(pts,1,1).^2,2)) ;
strainThFil = repmat( reshape( hd.InputData / raideurFil, [1,1,nbImgs] ),[3,1,1] ) ;
Comp.Strains =  seedHP.Strains -...
 repmat( ( seedHP.L0 - repmat( Comp.Lmortar, [3,1] ) ) ./ seedHP.L0 ,[1,1,nbImgs] ) .* strainThFil ;
Comp.meanStrain = reshape( meanNoNan(Comp.Strains,1), [ ],1 ) ;