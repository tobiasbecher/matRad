function matRad_plotLexicographicMixed(pln,ct,cst,resultGUIs)

for i = 1:numel(resultGUIs)
    for j = 1:numel(pln.originalPlans)
        figure
        phDose = pln.originalPlans(j).numOfFractions*resultGUIs{i}{j}.physicalDose;
        matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln.propStf(1).isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
        title(append(pln.originalPlans(j).radiationMode,int2str(i)))
    end
end

for i = 1:numel(resultGUIs)
    figure
    phDose = zeros(size(resultGUIs{1}{1}.physicalDose));
    for j = 1:numel(pln.originalPlans)
        phDose = phDose + pln.originalPlans(j).numOfFractions*resultGUIs{i}{j}.physicalDose;
    end
    matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln.propStf(1).isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
    title(append('Mixed',int2str(i)))
end

for i = 1:numel(resultGUIs)-1
    figure
    phDose = zeros(size(resultGUIs{1}{1}.physicalDose));
    phDose2 = zeros(size(resultGUIs{1}{1}.physicalDose));
    for j = 1:numel(pln.originalPlans)
        phDose = phDose + pln.originalPlans(j).numOfFractions*resultGUIs{i}{j}.physicalDose;
        phDose2 = phDose2 + pln.originalPlans(j).numOfFractions*resultGUIs{i+1}{j}.physicalDose;
    end
    phDose = phDose - phDose2;
    matRad_plotSliceWrapper(gca,ct,cst,1,phDose,3,round(pln.propStf(1).isoCenter(1,3)./ct.resolution.z),[],0.75,colorcube,[], [0 max([phDose(:)])],[]);
    title(append('Mixed',int2str(i)))
end