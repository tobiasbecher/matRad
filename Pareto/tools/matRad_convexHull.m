function [k,facets,normals,cs,dists,w] = matRad_convexHull(fVals,penalties)
% matRad function that takes the function values and associated penalties
% and returns the next penalty vector to run
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

%renormalize fVals so they are in same dimension for distance calculation
%in step 5

L = min(fVals,[],1);
U = max(fVals,[],1);
fVals = (fVals-L)./(U-L);

%calculate convex hull
[k,vol] = convhulln(fVals);
%%
%calculate normals for each facet(normal vectors are perpendicular to their respective facet)
%Done by solving P*n= 0 where P is a matrix with the vectors spanning the hyperplane and n the
%normal vector to be calculated

%initializing some objects that are returned by the function
normals = zeros(size(k));
cs  = zeros(size(k,1),1);
dists = zeros(size(k,1),1);
facets= zeros(size(k));
j = 0;
%% loop over all facets of convex hull
for i = 1:size(k,1)
    %% Step2: Calculate upper bounds from convex hull
    %vertices of facets
    ps = fVals(k(i,:),:);
    %calculate vectors spanning hyperplane through refPoint
    refPoint = ps(1,:);
    f = ps-refPoint; % vectors spanning hyperplan through ps(1,:) 
    zw = f(2:end,:);  % remove reference Point from 
    %
    normal = null(zw); %solve P*n=0 P Matrix
    
    %% check orientation of calculated normal vector
    %https://de.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra
    
    %get vertex not in facet and calculate orientation of facet
    idxs = (1:size(fVals,1));
    diffs = setdiff(idxs,k(i,:)); % find vertex not in facet
    
    %check if orientation and normal vector face in same direction
    orientationVector = fVals(diffs(1),:)-refPoint;
    hh = orientationVector*normal>0;
    orientation =(orientationVector*normal>0);    
    %flip orientation of normal vector if it goes in the opposite direction
    normal = normal*(2*orientation-1);
    
    
    %reject facet if normal has all negative components
    if all(normal<0)
        continue
    end
    j = j +1;
    %calculate constant for hyperplane equation (c = n*V)
    %n: normal Vector, V: Vector to point on hyperplane
    c = refPoint*normal;
    
    
    %% Step 3: Calculate lower bounds and Lower distal point (LDP)
    ws = penalties(k(i,:),:);
    %Upper bound hyperplane equation: ws(F-P)=0 find F for all 
    rs = sum(ws.*ps,2);
    %TODO: what happens if intersection in more than one point/no
    %intersection?
    LDP= linsolve(ws,rs);
    'a'
    %check that LDP is in bounding box of of PS, if not ignore facet
    if all(LDP<=ones(size(LDP))) && all(LDP>=zeros(size(LDP))) && (((refPoint-LDP')*normal) > 0)
        %% Step 5 Distance of LDP to upper bound
        dist = abs(LDP'*normal-c); % c: constant of upper bound
        facets(i,:) = k(i,:);
        cs(i,1) = c;
        normals(i,:) = normal;
        dists(i,1) = dist; 
    end
end
%% Step 6: Find new facet to run
%[argval,idx] = max(dists);
[A,I] = sort(dists,'descend');
%% Step 7: Find new weight vector
w = zeros(1,size(penalties,2));
%if all are nonnegative choose normal vector of facet as next penalties to
%run
found = false;
i = 1;
while ~found && i <= numel(I)
    idx = I(i);
    all(normals(idx,:)>=0)
    if all(normals(idx,:)>=0)  && ~all(normals(idx,:) == 0) && ~any(ismember(penalties,normals(idx,:),'rows'))
        w = normals(idx,:);
        found = true;
    elseif (~all(normals(idx,:) == 0))&& ~any(ismember(penalties,normals(idx,:),'rows')) %use maximally different vector
        PsPens = penalties(k(idx,:),:);
        %w = mean(PsPens);
        w = matRad_maxminVector(PsPens,mean(PsPens));
        found = true;
    end
    i = i+1;
end
facets
matRad_convexHullPlotHelper(fVals,facets,normals,penalties,w);
    
