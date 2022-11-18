function matRad_plotNormals(fVals,weights,VOIObj)
% matRad function that if orientation of penalty vectors align with PS
% input
%   fVals:      Array storing the function values
%   penalties:  Array storing the penalties used
% output
%   
%
% References
%   DOI:        10.1118/1.2335486
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
figure
L = min(fVals,[],1);
U = max(fVals,[],1);
fVals = (fVals-L)./(U-L);
%matRad_plotParetoSurface(fVals,weights,VOIObj);
plot(fVals(:,1),fVals(:,2),'Color','g','LineWidth',2,DisplayName = 'ParetoSurface');
hold on
plot(fVals(:,1),fVals(:,2),'*','Color','black','HandleVisibility', 'off');
%% Plot points 
size(fVals,1)
for i = 1:size(fVals,1)
    pts = [fVals(i,:);fVals(i,:)+weights(i,:)];
    wPlane = fliplr(weights);
    wPlane(:,1) = wPlane(:,1)*(-1);
    ptwPlane = [fVals(i,:);fVals(i,:)+wPlane(i,:);fVals(i,:)-wPlane(i,:),];
    if i==1
        %plot(pts(:,1),pts(:,2),'LineStyle','-.','Color','b', DisplayName = 'Penalty vectors')
        plot(ptwPlane(:,1),ptwPlane(:,2),'LineStyle','-','Color',[0,0,1,0.2],DisplayName = 'Lower Boundary')
    else
        %plot(pts(:,1),pts(:,2),'LineStyle','-.','Color','b', 'HandleVisibility', 'off')
        plot(ptwPlane(:,1),ptwPlane(:,2),'LineStyle','-','Color',[0,0,1,0.2],'HandleVisibility', 'off')
    end
end
legend()