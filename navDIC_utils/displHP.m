function hd = displHP( hd ) 

if ~isfield(hd,'Seeds') || isempty(hd.Seeds) ; return; end
nbIm = size(hd.Seeds(1).Displacements,3) ;
listSeed = {hd.Seeds(:).Name} ;
IDs = listdlg('PromptString','Select the Camera :',...
   'SelectionMode','Multiple',...
   'initialValue',1,...
   'ListString',listSeed) ;

seeds = hd.Seeds(IDs) ;
CamIDs = [seeds.CamIDs] ; 

% Choix du point de reference
fig = figure ;
nbPts = zeros(1,2) ; 
for i = 1:2
    nbPts(i) = size( seeds(i).Displacements, 1 ) ;
    ax = subplot(2,2,i) ;
    legtag = {} ;
    ind = {} ;
    for j = 1:nbPts
        plot( reshape( sum( seeds(i).Displacements(j,:,:).^2, 2 ).^.5, [],1 ) ) ;
        legtag{end+1} = [ 'Points numero ', num2str(j) ] ;
        ind{end+1} = num2str(j) ;
    end
    legend(legtag) ;
    ax = subplot(2,2,i+2) ;
    ax.YDir = 'reverse' ;
    ax.DataAspectRatio = [1,1,1] ;
    if iscell(hd.Images{1}{CamIDs(i)})
        image( repmat( hd.Images{1}{CamIDs(i)}{1}', [1,1,3] ) )
    else
        images( repmat( hd.Images{1}{CamIDs(i)} ), [1,1,3] )
    end

    plot( seeds(i).Points( :, 2 ), seeds(i).Points( :, 1 ),'x' )
    text(seeds(i).Points( :, 2 ), seeds(i).Points( :,1 ), ind )
end
%input( 'Press enter to continue :' ) ;

ind1 = listdlg('PromptString','Select point to calculate déplacement HP camera1:',...
                           'SelectionMode','multiple','ListSize',[250,250],...
                           'ListString',legtag);
ind2 = listdlg('PromptString','Select point to calculate déplacement HP camera2:',...
                           'SelectionMode','multiple','ListSize',[250,250],...
                           'ListString',legtag);

close(fig) ;

seeds(1).IDrefHP = ind1 ;
seeds(2).IDrefHP = ind2 ;

% calcul du mouvement hors plan 
if ~isfield( hd.Cameras(1), 'Properties') || isempty( hd.Cameras(1).Properties ) ; hd = setCameraProperties(hd,1) ; end
if ~isfield( hd.Cameras(2), 'Properties') || isempty( hd.Cameras(2).Properties ) ; hd = setCameraProperties(hd,2) ; end
X1 = hd.Cameras(CamIDs(1)).Properties.X ;
X2 = hd.Cameras(CamIDs(2)).Properties.X ;
Y1 = hd.Cameras(CamIDs(1)).Properties.Y ;
Y2 = hd.Cameras(CamIDs(2)).Properties.Y ;
Z1 = hd.Cameras(CamIDs(1)).Properties.Z ;
Z2 = hd.Cameras(CamIDs(2)).Properties.Z ;
P1 = hd.Cameras(CamIDs(1)).Properties.P ;
P2 = hd.Cameras(CamIDs(2)).Properties.P ;
pix2mm1 = 1/hd.Cameras(CamIDs(1)).Properties.pixObjratio ;
pix2mm2 = 1/hd.Cameras(CamIDs(2)).Properties.pixObjratio ;
f1 = hd.Cameras(CamIDs(1)).Properties.fpix(1) ;
f2 = hd.Cameras(CamIDs(2)).Properties.fpix(1) ;
do1 = hd.Cameras(CamIDs(1)).Properties.do ;
do2 = hd.Cameras(CamIDs(2)).Properties.do ;

U(:,:,1) = permute( meanNoNan( seeds(1).Displacements( ind1 ,: ,: ), 1 ), [3 2 1] ) ;
U(:,:,2) = permute( meanNoNan( seeds(2).Displacements( ind2 ,: ,: ), 1 ), [3 2 1] ) ;


ur1 = U(:,:,1) * [dot(Z2,X1); dot(Z2,Y1)] ;
ur2 = U(:,:,2) * [dot(Z1,X2); dot(Z1,Y2)] ;

h2 = zeros([nbIm,1]) ;
for i=1:4
    h1 = ur2 / f2 .* ( repmat(do2,[nbIm,1]) + h2 ) ; 
    h2 = ur1 / f1 .* ( repmat(do1,[nbIm,1]) + h1 ) ; 
end

seeds(1).Uhp = h1 ; %  U(:,:,2) * pix2mm2 * [dot(Z1,X2); dot(Z1,Y2)] ;
seeds(2).Uhp = h2 ; % U(:,:,1) * pix2mm1 * [dot(Z2,X1); dot(Z2,Y1)] ;

% Calcul de theta HP
% vectCam = zeros( [ nbIm, 3, 2 ] ) ;
% for i = 1:2
%     Pts = ones( [ nbPts(i), 3, nbIm ] ) * hd.Cameras( CamIDs( i ) ).Properties.fpix( 1 ) ;
%     Pts( :, 1:2, : ) = seeds(i).MovingPoints ;
%     vect = cross( Pts(2:2:end,:,:), Pts(1:2:end,:,:), 2) ; 
%     vect = vect./ repmat( sum( vect.^2, 2 ).^.5, [ 1, 3, 1] ) ;
%     % ranger les vecteurs avec directions principales(x ou y) positives
%     if abs( vect(1,1,1) ) < abs( vect(1,2,1) ) 
%         vect( repmat( vect(:,2,:) < 0, [ 1, 3, 1 ] ) ) = -vect( repmat( vect(:,2,:) < 0, [ 1, 3, 1 ] ) ) ;
%     else
%         vect( repmat( vect(:,1,:) < 0, [ 1, 3, 1 ] ) ) = -vect( repmat( vect(:,1,:) < 0, [ 1, 3, 1 ] ) ) ;
%     end
%     vectMean = permute(sumNoNan( vect, 1 ), [ 3, 2, 1 ] ) ;
%     vectCam(:,:,i) = vectMean./ repmat( sum( vectMean.^2, 2 ).^.5, [ 1, 3] ) ;
% end
% t3dPC(:,:,3) = cross( (P2\vectCam(:,:,2)')', (P1 \ vectCam(:,:,1)')', 2 ) ; 
% % if vect3D( 1, :, 3 ) * [ 0 ;0 ;1 ] < 0 
% %     vect3D(:,:,3) = -vect3D(:,:,3) ;
% % end
% t3dPC(:,:,1) = ( P1 * t3dPC(:,:,3)' )' ;
% t3dPC(:,:,2) = ( P2 * t3dPC(:,:,3)' )' ;

PtsPix(:,:,:,1) = seeds(1).MovingPoints ;
PtsPix(:,:,:,2) = seeds(2).MovingPoints ;
Cam(1) = hd.Cameras(CamIDs(1)).Properties ; 
Cam(2) = hd.Cameras(CamIDs(2)).Properties ; 

[t3dPC,~] = dirVect2cam(PtsPix , Cam) ;
% ranger les vecteurs avec directions principales(x ou y) positives
for i=1:2
    if abs( t3dPC(1,1,i) ) < abs( t3dPC(1,2,i) ) 
        if t3dPC(1,2,i)<0
            t3dPC(:,:,i) = -t3dPC(:,:,i) ;
        end
    else
        if t3dPC(1,1,i)<0
            t3dPC(:,:,i) = -t3dPC(:,:,i) ;
        end
    end
end

% Projection de l'angle 3D sur le plan et calcul de l'angle HP

vectProj = zeros( nbIm, 3, 2 ) ;

for i = 1:2
    if abs( t3dPC(1,2,i) ) > abs( t3dPC(1,1,i) )
        disp(['Camera ', num2str(i), ' : rotation autour de X']) ;
        vectProj(:,2:3,i) = t3dPC(:,2:3,i) ;
        
    elseif abs( t3dPC(1,2,i) ) < abs( t3dPC(1,1,i) )
        disp(['Camera ', num2str(i), ' : rotation autour de Y']) ;
        vectProj(:,1:2:3,i) = t3dPC(:,1:2:3,i) ;
    end
    vectProj(:,:,i) = vectProj(:,:,i) ./ repmat( sum( vectProj(:,:,i).^2, 2 ).^.5, [1,3] ) ;
    seeds(i).theta = asin( vectProj( :, 3, i ) ) ;
end

hd.Seeds(IDs) = seeds ;
