function [k,facets,normals,cs,dists,w] = matRad_convexHull(fVals,weights)

%renormalize fVals
L = min(fVals,[],1);
U = max(fVals,[],1);
fVals = (fVals-L)./(U-L);
%calculate convex hull
[k,vol] = convhulln(fVals);
%calculate normals for each facet Pn=0 (normal vector has to perpendicular to facet)
%Done by solving P*n= 0 where P is a matrix with the points as row vectors and n the
%vector to be calculated
normals = zeros(size(k));
facets= zeros(size(k));
cs  = zeros(size(k,1),1);
dists = zeros(size(k,1),1);
%% loop over all facets of convex hull
for i = 1:size(k,1)
    %% Step2: Calculate upper bounds
    %vertices of facets
    ps = fVals(k(i,:),:);
    refPoint = ps(1,:);
    f = ps-refPoint; % vectors spanning hyperplan through ps(1,:) 
    %f
    zw = f(2:end,:); 
    %
    normal = null(zw); %solve P*n=0 P Matrix
    
    %check orientation of calculated normal vector
    %https://de.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra
    
    %get vertex not in facet and calculate orientation of facet
    idxs = (1:size(fVals,1));
    diffs = setdiff(idxs,k(i,:));
    %check if orientation and normal vector face in same direction
    
    orientationVector = fVals(diffs(1),:)-refPoint;
    
    orientation =(orientationVector*normal>0);     
    %flip orientation of normal vector if it goes in the opposite direction
    %normal
    normal = normal*(2*orientation-1);
    %normal
    %reject facet if normal has all negative components
    if ~all(orientation)
        continue
    end
    %calculate constant for hyperplane equation (c = n*V)
    %n: normal Vector, V: Vector to point on hyperplane
    c = refPoint*normal;
    refPoint
    normal
    c
    cs(i,1) = c;
    
    normals(i,:) = normal;
    facets(i,:) = k(i,:);
    %% Step 3: Calculate lower bounds and Lower distal point (LDP)
    ws = weights(k(i,:),:);
    %Upper bound hyperplane equation: ws(F-P)=0 find F for all 
    rs = sum(ws.*ps,2);
    %TODO: what happens if intersection in more than one point/no
    %intersection?
    LDP= linsolve(ws,rs);
    %% Step 5 Distance of LDP to upper bound
    %LDP
    %normal
    %LDP'*normal;
    dist = abs(LDP'*normal-c);
    dists(i,1) = dist;
end
%% TO BE REMOVED 
matRad_convexHullPlotHelper(fVals,facets,normals,weights,cs);
[facets,dists]
%% Step 6: Find new facet to run
[argval,argmax] = max(dists);
%% Step 7: Find new weight vector
w = zeros(1,size(weights,2));
%if all are nonnegative choose normal vector of facet as next weights to
%run
normals
normals(argmax,:)
if all(normals(argmax,:)>=0)
    w = normals(argmax,:)
else 
    'TO BE IMPLEMENTED'
end
    
