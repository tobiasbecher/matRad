function mixedDose =  matRad_calcMixedDose(resultGUI,pln)
    mixedDose = zeros(size(resultGUI{1}.physicalDose)); 
    for i = 1:numel(pln.originalPlans)
        mixedDose = mixedDose + resultGUI{i}.physicalDose*pln.originalPlans(i).numOfFractions;
    end
end