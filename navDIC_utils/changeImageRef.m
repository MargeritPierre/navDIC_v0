function hd2 = changeImageRef( hd,n ) 
if ~isfield(hd,'Seeds') || isempty(hd.Seeds) ; return; end
nbSeeds = length( hd.Seeds ) ;
hd2 = hd ;
seeds = hd2.Seeds ;
nbIm = size( seeds(1).Displacements, 3 ) ;
for i = 1:nbSeeds
    seeds(i).Points = seeds(i).MovingPoints( :,:,n ) ; 
    seeds(i).MovingPoints = seeds(i).MovingPoints - repmat( seeds(i).MovingPoints(:,:,1) - seeds(i).MovingPoints( :,:,n ) ,[1,1, nbIm] ) ; 
    seeds(i).Displacements = seeds(i).MovingPoints - repmat( seeds(i).MovingPoints( :,:,n ), [ 1, 1, nbIm ] ) ;
    seeds(i).L0 = sum( ( seeds(i).MovingPoints( 2:2:end,:,n) - seeds(i).MovingPoints( 1:2:end,:,n) ).^2, 2 ).^.5 ;
    L = sum( ( seeds(i).MovingPoints(2:2:end,:,: ) - seeds(i).MovingPoints(1:2:end,:, :) ).^2, 2 ).^.5 ;
    seeds(i).Strains  = ( L - repmat( seeds(i).L0, [ 1, 1, nbIm ] ) ) ./ repmat( seeds(i).L0, [ 1, 1, nbIm ] ) ;
end

list = {'oui','non'} ;
indx = listdlg('PromptString','Voulez-vous supprimer les images anterieures? ',...
                           'SelectionMode','single','ListSize',[250,50],...
                           'ListString',list);
if indx == 1
    for i = 1:nbSeeds
        seeds(i).Strains = seeds(i).Strains(:,:,n:end) ;
        seeds(i).Displacements = seeds(i).Displacements(:,:,n:end) ;
        seeds(i).MovingPoints = seeds(i).MovingPoints(:,:,n:end) ;
        hd2.InputData = hd.InputData(n:end,end) ;
    end
end

hd2.Seeds = seeds ;
end