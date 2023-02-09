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

%STEP 1: Done outside of function

%Preprocessing step

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
rejected_facets_normal = zeros(size(k));
rejected_facets_upper = zeros(size(k));
rejected_facets_lower = zeros(size(k));
rejected_facets_dist = zeros(size(k));
j = 0;
jj = 0;
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
    normal = null(zw); %solve P*n=0 P Matrix for calculation of normal vectors
    
    %% check orientation of calculated normal vector and check orientation
    %https://de.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra
    
    %get vertex not in facet and calculate orientation offacet
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
        rejected_facets_normal(i,:) = k(i,:);
        continue
    end
    j = j +1;
    %calculate constant for hyperplane equation (c = n*V)
    %n: normal Vector, V: Vector to point on hyperplane
    c = refPoint*normal; %(used for Step 4)
    
    
    %% Step 3: Calculate lower bounds and Lower distal point (LDP)
    ws = penalties(k(i,:),:);
    %Upper bound hyperplane equation: ws(F-P)=0 find F for all 
    rs = sum(ws.*ps,2);
    %TODO: what happens if intersection in more than one point/no
    %intersection?
    LDP= linsolve(ws,rs);
    'a'
    
    %check that LDP is in bounding box of of PS, if not ignore facet
    if ~all(LDP<=ones(size(LDP))*10)
        rejected_facets_upper(i,:) = k(i,:);
    elseif ~all(LDP>=zeros(size(LDP)))
        rejected_facets_lower(i,:) = k(i,:);
    elseif ~((refPoint-LDP')*normal) > 0
        rejected_facets_dist(i,:) = k(i,:);
    end 
    if all(LDP>=zeros(size(LDP))) && (((refPoint-LDP')*normal) > 0) && all(LDP<=ones(size(LDP))*10) 
        %% Step 4 Distance of LDP to upper bound
        dist = abs(LDP'*normal-c) ; % c: constant of upper bound
        facets(i,:) = k(i,:);
        cs(i,1) = c;
        normals(i,:) = normal;
        dists(i,1) = dist; 
        jj = jj +1;
    else
        rejected_facets(i,:) = k(i,:);
    end
end
jj/size(k,1)
%% Step 5 and 6: Find new facet to run
%[argval,idx] = max(dists);
[A,I] = sort(dists,'descend');
w = zeros(1,size(penalties,2));
%if all are nonnegative choose normal vector of facet as next penalties to
%run
found = false;
i = 1;
while ~found && i <= numel(I)
    accuracy = 2;
    idx = I(i);
    norm = normals(idx,:);
    PsPens = penalties(k(idx,:),:);
    %norm = norm/sum(norm);
    all(normals(idx,:)>=0)
    %if all(norm>=0)  && ~all(norm == 0) && ~any(ismember(round(penalties,accuracy),round(norm,accuracy),'rows'))
    %    w = normals(idx,:);
    %    zwvec = w/sum(w);
    %    found = true;
    if (~all(normals(idx,:) == 0))&& ~any(ismember(penalties,normals(idx,:),'rows')) %use maximally different vector
        w = mean(PsPens);
       
        %maxminw0 = rand(size(mean(PsPens)));
        %maxminw0 = maxminw0/sum(maxminw0);
        %w = matRad_maxminVector(PsPens,maxminw0);
        %w = matRad_maxminVector(PsPens,w);
        %zwvec = w;
        w = w/sqrt((sum(w.^2)));
        if all(w>=0)  && ~all(w == 0) && ~any(ismember(round(penalties,accuracy),round(w,accuracy),'rows'))
            found = true;
        end
    end
    i = i+1;
end
w = w/sqrt((sum(w.^2)))%renormalize

%matRad_convexHullPlotHelper(fVals,facets,normals,penalties,w);
%%
%{
figure
fill3(PsPens(:,1),PsPens(:,2),PsPens(:,3),'red')
hold on 
scatter3(zwvec(1),zwvec(2),zwvec(3))
%}
%%
figure
trisurf(k,fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','cyan')
hold on

trisurf(rejected_facets_normal(~all(rejected_facets_normal==0,2),:),fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','red')
trisurf(rejected_facets_upper(~all(rejected_facets_upper==0,2),:),fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','yellow')
trisurf(rejected_facets_upper(~all(rejected_facets_upper==0,2),:),fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','red')
trisurf(rejected_facets_lower(~all(rejected_facets_lower==0,2),:),fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','black')
trisurf(rejected_facets_dist(~all(rejected_facets_dist==0,2),:),fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','blue')
%%
facetPoints = fVals(k(idx,:),:)
%%
fill3(facetPoints(:,1),facetPoints(:,2),facetPoints(:,3),'green')
%%

zlim([0,1])   
xlim([0,1])
ylim([0,1])
hold on
xlabel('x')
ylabel('y')
zlabel('z')
'A'
%hold on
%point = fVals(4,:);
%normalPlot = penalties(4,:);

%# a plane is a*x+b*y+c*z+d=0
%# [a,b,c] is the normal. Thus, we have to calculate
%# d and we're set
%d = -point*normalPlot'; %'# dot product for less typing

%# create x,y
%[xx,yy]=ndgrid(0:1,0:1);

%# calculate corresponding z
%z = (-normalPlot(1)*xx - normalPlot(2)*yy - d)/normalPlot(3);

%# plot the surface
%surf(xx,yy,z)