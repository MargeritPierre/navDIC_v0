function [strainNerror, strainSerror] = systemError(StrainVect,StrainReal)

% Fonction qui calcule le bruit de mesure à partir de plusieurs images
% prises successivement en deplaçant la camera avec la platine micrometrique.


E = StrainVect;

E(isnan(E)) = 1i ;
strainNmean = sum(real(E),1)./sum(imag(E)==0,1) ;

Emean = repmat(strainNmean, [length(E) 1]);
Emean(real(E)==0) = 1i;

strainNerror = sum(abs(real(E-Emean)),1)./sum(imag(E)==0,1);

strainSerror = abs(strainNmean - StrainReal);

disp(['l erreur de bruit en Exx est de ',num2str(strainNerror(1,1))]);
disp(['l erreur de bruit en Eyy est de ',num2str(strainNerror(1,2))]);
disp(['l erreur de bruit en Exy est de ',num2str(strainNerror(1,3))]);

disp(['l erreur systematique en Exx est de ',num2str(strainSerror(1,1))]);
disp(['l erreur systematique en Eyy est de ',num2str(strainSerror(1,2))]);
disp(['l erreur systematique en Exy est de ',num2str(strainSerror(1,3))]);

end