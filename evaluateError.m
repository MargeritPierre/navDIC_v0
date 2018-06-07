[strainNerror1, strainSerror1] = systemError(hd.Seeds.Strains(:,:,5),10^(-5));

[displNerror2, strainNerror2] = noiseError(hd.Seeds.Displacements(:,:,5), hd.Seeds.Strains(:,:,5));

mean(hd.Seeds.Displacements(:,:,5),1)

mean(hd.Seeds.Strains(:,:,5),1)

strainNerror = zeros(hd.CurrentFrame, 3);
strainSerror = zeros(hd.CurrentFrame, 3);

for i=1:hd.CurrentFrame
    [strainNerror(i,:), strainSerror(i,:)] = systemError(hd.Seeds.Strains(:,:,i),10^(-5));
    
end

realStrain = 10^(-5);

save([hd.WorkDir.Path, 'StrainError.mat'], 'strainNerror', 'strainSerror','realStrain');