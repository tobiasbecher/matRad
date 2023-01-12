function wmin = matRad_maxminVector(pens,w0)
    % matRad function that creates a penalty vector maximally different to other
    % vectors of a facet
    %
    % input
    %   pens:       Matrix storing the penalty values of the facet
    %   w0:         Initial point (could be initialize
    %
    % output
    %   wmin:       Maximally different vector
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2015 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lb = zeros(1,size(pens,1));
    ub = ones(1,size(pens,1));
    %options = optimoptions('fmincon',...
        %'OptimalityTolerance');
        %'ConstraintTolerance',1);
    %ones(1,size(pens,1))
    wmin = fminimax(@(al) matRad_allVectorDiff(al,pens),...
                    w0,... % Starting Point
                    [],[],...% Linear Constraints: sum(al)=1;
                    ones(1,size(pens,1)),1,... % Also no linear inequality constraints
                    lb,ub);%,options);
    wmin
    pens
    matRad_vectorDiff(wmin,pens)
    wmin = wmin*pens;