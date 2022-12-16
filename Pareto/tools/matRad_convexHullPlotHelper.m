function returnStruct = matRad_convexHullPlotHelper(fVals,facets,normals,weights,nw)
% matRad inverse pareto planning wrapper function
% 
% call
% output
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
[fVals,weights];
[k,vol] = convhulln(fVals);
figure;
plot(fVals(:,1),fVals(:,2),"diamond",color = "k");
hold on 
facets
facets = facets(any(facets,2),:);
facets
normals = normals(any(normals,2),:);
difs = fVals(facets(:,1),:)-fVals(facets(:,2),:);
medians = fVals(facets(:,1),:)-1/2*difs;
%%
plot(medians(:,1),medians(:,2),'*','HandleVisibility', 'off');
plot(fVals(:,1),fVals(:,2),'.','HandleVisibility', 'off');
fVals
fVals(facets',1)
fVals(facets',2)
%plot(fVals(k',1),fVals(k',2),'color','k',DisplayName = 'Convex Hull',LineWidth = 1.3);
plot(fVals(facets',1),fVals(facets',2),'color','k',DisplayName = 'Considered facets',LineWidth = 1.3);

%plot normals of convex hulls

for i = 1:size(facets,1)
    vec = [medians(i,:);medians(i,:)+normals(i,:)];
    if i==1
        plot(vec(:,1),vec(:,2),'LineStyle',':','Color','r',DisplayName = 'Facet normals');
    else
        plot(vec(:,1),vec(:,2),'LineStyle',':','Color','r','HandleVisibility', 'off');
    end
end
%%

nidx = find(ismember(normals, nw,'rows'));
vec = [medians(nidx,:);medians(nidx,:)+normals(nidx,:)];
plot(vec(:,1),vec(:,2),'LineStyle','-','Color','r',DisplayName = 'Next normal',LineWidth = 2);

%%
%

for i = 1:size(fVals,1)
    pts = [fVals(i,:);fVals(i,:)+weights(i,:)];
    
    wPlane = fliplr(weights);
    wPlane(:,1) = wPlane(:,1)*(-1);
    ptwPlane = [fVals(i,:);fVals(i,:)+wPlane(i,:);fVals(i,:)-wPlane(i,:),];
    if i==1
        plot(pts(:,1),pts(:,2),'LineStyle','-.','Color','b', DisplayName = 'Penalty vectors')
        plot(ptwPlane(:,1),ptwPlane(:,2),'LineStyle','--','Color','g',DisplayName = 'Lower Boundary')
    else
        plot(pts(:,1),pts(:,2),'LineStyle','-.','Color','b', 'HandleVisibility', 'off')
        plot(ptwPlane(:,1),ptwPlane(:,2),'LineStyle','--','Color','g','HandleVisibility', 'off')
    end
end
%%


legend()
%}
