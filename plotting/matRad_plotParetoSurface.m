function matRad_plotParetoSurface(fInd,penGrid)
% matRad function that plots a colour coded Pareto Surface. 
%
% call
%   cBarHandle = matRad_plotColorbar(axesHandle,cMap,window,varargin)
%
% input
%   axesHandle  handle to axes the colorbar will be attached to
%   cMap        corresponding colormap
%   window      colormap window (corresponds to clim)
%   varargin    additional key-value pairs that will be forwarded to the
%               MATLAB colorbar(__) call
%
% output
%   cBarHandle  handle of the colorbar object
%
% References
%   -
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
fig2 = figure;
switch size(fInd,2)
    case 2
        penGrid = [penGrid,zeros(size(penGrid,1),1)];
        scatter(fInd(:,1),fInd(:,2),[],penGrid,'filled')
        xlabel('Objective function value 1');
        ylabel('Objective function value 2');
    case 3
        scatter3(fInd(:,1),fInd(:,2),fInd(:,3),[], penGrid,'filled')
        xlabel('Objective function value 1');
        ylabel('Objective function value 2');
        zlabel('Objective function value 3');
    otherwise
        warning(['Number of objectives for Pareto Analysis not suited for Plot!']);
end
