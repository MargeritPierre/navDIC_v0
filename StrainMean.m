function [strainMean] = StrainMean(StrainVect)

% Fonction qui calcule le bruit de mesure à partir de plusieurs images
% prises successivement en deplaçant la camera avec la platine micrometrique.


E = StrainVect;

E(isnan(E)) = 1i ;
strainMean = sum(real(E),1)./sum(imag(E)==0,1) ;

end