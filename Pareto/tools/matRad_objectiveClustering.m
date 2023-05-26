function [obj1,obj2] = matRad_objectiveClustering(fVals)
    % matRad function that clusters objectives based on their spearman rank
    % correlation
    %
    % input
    %  fVals:   All calculated objective function values
    % output
    %   obj1:   Index of objective 1 used for grouping
    %   obj2:   Index of objective 2 used for grouping
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


    rho = triu(corr(fVals,'Type','Spearman'),1); %Only use upper matrix f
    [M,I] = max(rho,[],'all');
    thresold = 0.5;
    if rho(I) > threshold;
       %convert index 
       [obj1,obj2] =  ind2sub(size(rho),I);



