function fValsMod = matRad_generateDummyPoints(fVals) 
    % matRad Function that generates Dummy Points for a set of points
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
    
theta = 0.1;
m = size(fVals,2);
n = size(fVals,1);

U = max(fVals,[],1)*m+theta;

fValsMod = zeros(size(fVals,1)*size(fVals,2),size(fVals,2));
for i = 1:n
    temp = repmat(fVals(i,:),m,1);
    for j = 1:numel(U)
        temp(j,j) = U(j);
    end
    fValsMod((i-1)*m+1:i*m,:) = temp;
end
fValsMod = [fVals;fValsMod];