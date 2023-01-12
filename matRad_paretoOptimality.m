function paretoOptimalPoints = matRad_paretoOptimality(fVals)
    paretoOptimalPoints = zeros(size(fVals));
    for i = 1:size(fVals,1)
        j = 1;
        pOptimal = true;
        while pOptimal && j<= size(fVals,1)
            %check if one point dominates the other
            if all(fVals(j,:) < fVals(i,:))
                pOptimal = false;
            end
            j = j+1;
        end
        if pOptimal
            paretoOptimalPoints(i,:) = fVals(i,:);
        end
    end
    %remove non zero lines (possible issues if "all zeros" is an actual solution but should never occur for radiotherapy results)
    paretoOptimalPoints = paretoOptimalPoints(all(paretoOptimalPoints>0,2),:) 