% importer points logiciels externes

for k = 1:length(hd.Cameras)
    clear pts
    [file,path] = uigetfile('.csv',['Selectionner le fichier de points de la camera ', hd.Cameras(k).Name,': ']) ;
    impPts = csvread([path,file],2,1) ;
    hd = setCameraProperties(hd,k) ;
    rezY = hd.Cameras(k).Properties.py ;
    ordre = input('Entrer 1 si les colonnes X et Y sont altern?es, 0 sinon : ') ;
    for i = 1:size(impPts, 1)
        for j = 1:size(impPts, 2)/2
                if ordre
                    pts(j,1,i) = impPts(i,(j-1)*2 + 1) ;
                    pts(j,2,i) = rezY - impPts(i,(j-1)*2+2 ) ;
                else
                    pts(j,1,i) = impPts(i,j) ;
                    pts(j,2,i) = rezY - impPts(i,size(impPts, 2)/2 + j ) ;
                end
        end
    end
    pts(pts==0)=NaN;
    PtsPix(:,:,:,k) = pts ;
end


%%
imRef = 1 ;
Cam(1) = hd.Cameras(1).Properties ;
Cam(2) = hd.Cameras(2).Properties ;
PtsPixComp = PtsPix(1:2,:,:,:);
PtsPixFil = PtsPix(3:4,:,:,:);

t3DComp = dirVect2cam(PtsPixComp , Cam) ;
t3DFil = dirVect2cam(PtsPixFil , Cam) ;

PtsPixCompCor = PixCor2cam(PtsPixComp, Cam, imRef) ;
PtsPixFilCor = PixCor2cam(PtsPixFil, Cam, imRef) ;
%PtsPixFil(1,:,1,1) 
%(PtsPixFilCor(1,1:2,1,1) + Cam(1).PPL')

%%

nbIm = size(PtsPixFil,3) ;
Lech = 10 ;

for i = 1:length(Cam)
% pixel apparent
    L0Comp(1,i) = sum( ( PtsPixComp(2,:,imRef,i) - PtsPixComp(1,:,imRef,i) ).^2, 2 ).^0.5 ;
    LComp(:,i) = sum( ( squeeze(PtsPixComp(2,:,:,i))' - squeeze(PtsPixComp(1,:,:,i))' ).^2, 2 ).^0.5 ;
    dLComp(:,i) = LComp(:,i) - repmat(L0Comp(1,i),[nbIm 1]) ;
    defComp(:,i) = dLComp(:,i) ./ repmat(L0Comp(1,i),[nbIm 1]) ;
    
    L0Fil(1,i) = sum( ( PtsPixFil(2,:,imRef,i) - PtsPixFil(1,:,imRef,i) ).^2, 2 ).^0.5 ;
    LFil(:,i) = sum( ( squeeze(PtsPixFil(2,:,:,i))' - squeeze(PtsPixFil(1,:,:,i))' ).^2, 2 ).^0.5 ;
    dLFil(:,i) = LFil(:,i) - repmat(L0Fil(1,i),[nbIm 1]) ;
    defFil(:,i) = dLFil(:,i) ./ repmat(L0Fil(1,i),[nbIm 1]) ;
% pixel corrig?s
    L0CompCor(1,i) = sum( ( PtsPixCompCor(2,:,imRef,i) - PtsPixCompCor(1,:,imRef,i) ).^2, 2 ).^0.5 ;
    LCompCor(:,i) = sum( ( squeeze(PtsPixCompCor(2,:,:,i))' - squeeze(PtsPixCompCor(1,:,:,i))' ).^2, 2 ).^0.5 ;
    dLCompCor(:,i) = LCompCor(:,i) - repmat(L0CompCor(1,i),[nbIm 1]) ;
    defCompCor(:,i) = dLCompCor(:,i) ./ repmat(L0CompCor(1,i),[nbIm 1]) ;
    
    L0FilCor(1,i) = sum( ( PtsPixFilCor(2,:,imRef,i) - PtsPixFilCor(1,:,imRef,i) ).^2, 2 ).^0.5 ;
    LFilCor(:,i) = sum( ( squeeze(PtsPixFilCor(2,:,:,i))' - squeeze(PtsPixFilCor(1,:,:,i))' ).^2, 2 ).^0.5 ;
    dLFilCor(:,i) = LFilCor(:,i) - repmat(L0FilCor(1,i),[nbIm 1]) ;
    defFilCor(:,i) = dLFilCor(:,i) ./ repmat(L0FilCor(1,i),[nbIm 1]) ;
    
    defCompCorReel(:,i) = defCompCor(:,i) * L0CompCor(1,i) / (Lech * Cam(i).pixObjratio) -...
        defFilCor(:,i) * ( L0CompCor(1,i) - (Lech * Cam(i).pixObjratio) ) / (Lech * Cam(i).pixObjratio) ;
end

%% Calcul des erreurs 
drot = [0.5/180*pi,0.5/180*pi] ; 
df = 0 ;
errDef = projError(PtsPixCompCor,Cam,drot,df) ;

%%

figure ;
plot(defFil)


%% Plots 
%path = uigetdir('Choisissez le fichier d''enregistrement des plots : ') ;
mkdir(hd.WorkDir.Path,'/output-matLab-GOM')
ruin = 200 ;
path = [ hd.WorkDir.Path,'/output-matLab-GOM','/'] ;

save([path,'defvalue'],'defFil','defComp','defFilCor','defCompCor','defCompCorReel')
% Trace defomation non corrig?e vs corrig?e cam 1 
%fil
fig(1) = figure ;
plot(defFil(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
plot(defFilCor(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{f} non corrigees','\epsilon_{f} corrigees'},'Interpreter','tex','Location','southeast') ;
title('Force vs Deformations du fil camera 1')

savefig(fig(1),[path,'Force vs Deformations du fil camera 1']) ;

fig(2) = figure ;
plot(defFil(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
plot(defFilCor(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{f} non corrigees','\epsilon_{f} corrigees'},'Interpreter','tex','Location','southeast') ;
title('Force vs Deformations du fil camera 2')

savefig(fig(2),[path,'Force vs Deformations du fil camera 2']) ;

% comp
fig(3) = figure ;
plot(defComp(imRef:ruin,1),hd.InputData(imRef:ruin)-hd.InputData(1))
plot(defCompCor(imRef:ruin,1),hd.InputData(imRef:ruin)-hd.InputData(1))
plot(defCompCorReel(imRef:ruin,1),hd.InputData(imRef:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{app} non corrigees','\epsilon_{app} corrigees','\epsilon_{c}'},...
    'Location','southeast','Interpreter','tex') ;
title('Force vs Deformations du composite camera 1')

savefig(fig(3),[path,'Force vs Deformations du composite camera 1']) ;

% comp
fig(4) = figure ;
plot(defComp(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
plot(defCompCor(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
plot(defCompCorReel(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{app} non corrigees','\epsilon_{app} corrigees','\epsilon_{c}'},...
    'Location','southeast','Interpreter','tex') ;
title('Force vs Deformations du composite camera 2')

savefig(fig(4),[path,'Force vs Deformations du composite camera 2']) ;


% 
fig(5) = figure ;
plot(defFilCor(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
plot(defCompCorReel(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
plot(defFilCor(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
plot(defCompCorReel(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{f} camera 1','\epsilon_{c} camera 1','\epsilon_{f} camera 2','\epsilon_{c} camera 2'},...
    'Location','southeast','Interpreter','tex') ;
title('Force vs Deformations corrig?es des deux cameras')

savefig(fig(5),[path,'Force vs Deformations corrig?es des deux cameras']) ;

% 
fig(6) = figure ;
plot(mean([defFilCor(1:ruin,1),defFilCor(1:ruin,2)],2),hd.InputData(1:ruin)-hd.InputData(1))
plot(mean([defCompCorReel(1:ruin,1),defCompCorReel(1:ruin,2)],2),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{f}','\epsilon_{c}'},...
    'Location','southeast','Interpreter','tex') ;
title('Force vs deformation corrigees, moyennes des deux cameras')

savefig(fig(6),[path,'Force vs deformation corrigees, moyennes des deux cameras']) ;

%%

k = 1 ;
for i = 1:length(hd.Cameras)
    ax(i) = subplot(1,2,i) ;
%      axes('Parent', fig, 'Units', 'normalized', 'PlotBoxAspectRatio' ,  [1 1 1], 'OuterPosition', pos(i,:), ...
%         'DataAspectRatio' , [1 1 1], 'XLim' , [0 , Cam(i).px] , 'YLim' , [0 , Cam(i).py] ) ;

    plotcam = plot( PtsPixFil(:,1,k,i), PtsPixFil(:,2,k,i),'*' ) ;%,'Parent',ax(i) ) ;
    plotcam = plot( PtsPixFilCor(:,1,k,i) + Cam(i).PPL(1) , PtsPixFilCor(:,2,k,i)+Cam(i).PPL(2),'*' ) ;
    axis equal
    xlim([0 , hd.Cameras(i).Properties.px]) ;
    ylim([ 0, hd.Cameras(i).Properties.py ]) ;
    ax(i).YDir='reverse' ;
    title(['Camera', num2str(i)] ) ;
    hold off
end




%% plot
k = 1 ;
ax = axes ;
plot(pts(:,1,k),pts(:,2,k),'*')
plot(pts(:,1,250),pts(:,2,250),'*')
axis equal
xlim([0 , hd.Cameras(2).Properties.px ]) ;
ylim([ 0, hd.Cameras(2).Properties.py  ]) ;
ax.YDir='reverse' ;
title('Imported points') ;
hold off

%%
fig = figure ;
%plot(defFil(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
%plot(defComp(1:ruin,1),hd.InputData(1:ruin)-hd.InputData(1))
plot(defFilCor(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
plot(defCompCorReel(1:ruin,2),hd.InputData(1:ruin)-hd.InputData(1))
xlabel('Deformation') ;
ylabel('Force (N)') ;
legend({'\epsilon_{f} camera 1','\epsilon_{c} camera 1','\epsilon_{f} camera 2','\epsilon_{c} camera 2'},...
    'Location','southeast','Interpreter','tex') ;
title('Force vs Deformations corrig?es des deux cameras')