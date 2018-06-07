function [displNerror, strainNerror] = noiseError(displVect, StrainVect)

% Fonction qui calcule le bruit de mesure à partir de plusieurs images
% prises successivement dun objet fixe.

U = displVect;

U(isnan(U)) = 1i ;
displNerror = sum(abs(real(U)),1)./sum(imag(U)==0,1) ;


disp(['l erreur de bruit en Ux est de ',num2str(displNerror(1)),'pixels']);
disp(['l erreur de bruit en Uy est de ',num2str(displNerror(2)),'pixels']);

E = StrainVect;

E(isnan(E)) = 1i ;
strainNerror = sum(abs(real(E)),1)./sum(imag(E)==0,1) ;

disp(['l erreur de bruit en Exx est de ',num2str(strainNerror(1))]);
disp(['l erreur de bruit en Eyy est de ',num2str(strainNerror(2))]);
disp(['l erreur de bruit en Exy est de ',num2str(strainNerror(3))]);

end