function resultGUI = matRad_calcCubesDoseGridMixedWrapper(resultGUIMixed,dij,scenNum)
    for mod = 1:numel(dij.original_Dijs)
        resultGUI{mod} = matRad_calcCubesDoseGrid(resultGUIMixed{mod}.w,dij.original_Dijs{mod},scenNum);
    end
    

end
