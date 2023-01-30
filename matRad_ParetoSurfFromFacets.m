function [k,facets] = matRad_ParetoSurfFromFacets(fVals)
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
    if any(normal<0)
        continue
    end
    facets(i,:) = k(i,:);
end
'a'