function f = matRad_objectiveFunction2(optiProb,w,dij,cst)
    % matRad IPOPT objective function wrapper
    %
    % call
    %   f = matRad_objectiveFuncWrapper(optiProb,w,dij,cst)
    %
    % input
    %   optiProb: matRad optimization problem
    %   w:        beamlet/ pencil beam weight vector
    %   dij:      matRad dose influence struct
    %   cst:      matRad cst struct
    %   scenario: index of dij scenario to consider (optional: default 1)
    %
    % output
    %   f: objective function value
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2016 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    fIndv = matRad_objectiveFunctions(optiProb,w,dij,cst);
    %fIndv should be an mxn matrix where m is the number of objectives and n should be the number of scenarios used

    switch optiProb.normalizationScheme.type
        case 'UL' %does not make sense together with robustness or scenarios? -> Add check?
            fIndv = optiProb.normalizeObjectives(fIndv);%ISSUES IF COWC is used -
    end

    %TODO: sum up over scenarios
    

    %COWC calculation

    %get indices of objectives and their respective penalties
    COWCidxs = [];
    idxs = [];
    penalties = [];
    COWCpenalties = [];
    f = 0;

    for i = 1:numel(optiProb.objectives)
        if strcmp(optiProb.objectives{i}.robustness,'COWC') 
            COWCidxs = [COWCidxs,i];
            COWCpenalties = [COWCpenalties,optiProb.objectives{i}.penalty];
        else
            idxs = [idxs,i];
            penalties = [penalties,optiProb.objectives{i}.penalty];
        end
    end
    if isempty(COWCidxs)
        f_COWC  = [0];
    else
        f_COWCs = fIndv(COWCidxs,:);
        f_COWC = COWCpenalties*f_COWCs; 
    end%sum over scenarios respect penalties
        
    fMax = max(f_COWC(:)); 
    if fMax > 0
        switch optiProb.useMaxApprox
            case 'logsumexp'
                fMax = optiProb.logSumExp(f_COWC);
            case 'pnorm'
                fMax = optiProb.pNorm(f_COWC,numel(useScen));
            case 'none'
                fMax = max(f_COWC);
            case 'otherwise'
                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                fMax = max(f_COWC);
        end
    end
    
    f = f + fMax;

    %now sum up over other robustness scenarios
    if isempty(idxs)
        f_wo_COWC = [0];
    else
        f_wo_COWCs = fIndv(idxs,:);
        f_wo_COWC = sum(penalties*f_wo_COWCs,2); %sum over scenarios?
    end
    f = f + f_wo_COWC;
end