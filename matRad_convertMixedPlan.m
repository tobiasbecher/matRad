function resultGUIVis = matRad_convertMixedPlan(resultGUI,pln,bound)
    %for now only physicalDose
    if nargin < 3
        bound = 1e9;
    end
    resultGUIVis.physicalDose = zeros(size(resultGUI{1}.physicalDose)); 
    for i = 1:numel(pln.originalPlans)
        phDose = resultGUI{i}.physicalDose*pln.originalPlans(i).numOfFractions;
        phDose(phDose> bound) = bound;
        resultGUIVis.physicalDose = resultGUIVis.physicalDose + phDose;
        resultGUIVis.(['physicalDose', pln.originalPlans(i).radiationMode]) = phDose;
    end

    %calculate relative values
    if numel(resultGUI) > 2
        for i = 1:numel(pln.originalPlans)
            phDose = resultGUI{i}.physicalDose*pln.originalPlans(i).numOfFractions;
            phDose(phDose> bound) = bound;
            resultGUIVis.(['relativeDose', pln.originalPlans(i).radiationMode]) = phDose./resultGUIVis.physicalDose;
        end
    
        resultGUIVis.FracDose =  resultGUI{2}.physicalDose*pln.originalPlans(2).numOfFractions;
        resultGUIVis.FracDose =  resultGUIVis.FracDose./ resultGUI{1}.physicalDose*pln.originalPlans(1).numOfFractions;
        resultGUIVis.FracDose(resultGUIVis.FracDose >1e9) = 1e9;
    end
end