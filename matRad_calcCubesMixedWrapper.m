function resultGUI = matRad_calcCubesMixedWrapper(resultGUIMixed,dij,scenNum,rescaleVar)
    for mod = 1:numel(dij.original_Dijs)
        resultGUI{mod} = matRad_calcCubes(resultGUIMixed{mod}.w,dij.original_Dijs{mod},scenNum,rescaleVar);
    end
    

end
