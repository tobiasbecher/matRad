function [cl,cu] = matRad_getConstraintBounds2(optiProb,cst)
    % matRad IPOPT get constraint bounds wrapper function
    % 
    % call
    %   [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
    %
    % input
    %   cst:            matRad cst struct
    %
    % output
    %   cl: lower bounds on constraints
    %   cu: lower bounds on constraints
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
    
    
    % Initialize bounds
    cl = [];
    cu = [];
    
    % compute objective function for every VOI.
    for i = 1:size(optiProb.constridx,1)	
        obj = optiProb.constraints{i};
        curConidx = optiProb.constridx(i,1);

        cl = [cl;obj.lowerBounds(numel(cst{curConidx,4}{1}))];
        cu = [cu;obj.upperBounds(numel(cst{curConidx,4}{1}))];
            
    end

       
    
    