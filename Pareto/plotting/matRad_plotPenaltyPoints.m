function matRad_plotPenaltyPoints(penPoints)
% matRad function to plot the penalty Vectors - works only for 2 or 3
% dimensional objectives
%
%
% input
%   penPoints    Matrix storing the penalty vectors
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
switch size(penPoints,2)
    case 2
        scatter(penPoints(:,1),penPoints(:,2),'filled','MarkerFaceColor',[0.4660 0.6740 0.1880])
        xlabel('Penalty 1');
        ylabel('Penalty 2');
    case 3
        scatter3(penPoints(:,1),penPoints(:,2),penPoints(:,3),[], penPoints(:,3),'filled')
        colormap(gca,"summer")
        xlabel('Penalty 1');
        ylabel('Penalty 2');
        zlabel('Penalty 3');
    case 4
        scatter3(penPoints(:,1),penPoints(:,2),penPoints(:,3),[], penPoints(:,4),'filled')
        h = colorbar()
        colormap(gca,"summer")
        xlabel('Penalty 1');
        ylabel('Penalty 2');
        zlabel('Penalty 3');
        set(get(h,'title'),'string','Penalty 4');
    otherwise
        warning(['Number of objectives for Pareto Analysis not suited for Plot!']);
end
