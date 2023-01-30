function [k,facets,normals,distfacet,w] =matRad_PGEN(fVals,penalties,distVal)
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

%normalize
L = min(fVals,[],1);
U = max(fVals,[],1);
fVals = (fVals-L)./(U-L);



%calculate convexHull
[k,vol] = convhulln(fVals);



%initializing

dists = zeros(size(k,1),1);
facets= zeros(size(k));
%we need the normals of the facets
normals = zeros(size(k));

%loop over all facets in convex hull
for i =1:size(k,1)
    
    facetPoints = fVals(k(i,:),:); %points that span the given facet
    %choose a reference point of facet that hyperplane is "build on"
    refPoint = facetPoints(1,:);
    hyperplaneVectors = facetPoints-refPoint; % calculate difference of reference point to all other points on facet
    spanningVectors = hyperplaneVectors(2:end,:); %reference point not needed for calculation of normal
    
    %calculate the normal of the hyperplane by solving V*n = 0 where V is
    %the matrix with the vectors spanning the hyperplane
    normal = null(spanningVectors);

    %we want to check if the components of the normal contain negative
    %components. For that one has to check the orientation first
    
    %choose a point of the PS that is not part of the hyperplane currently
    %investigated
    idxs = (1:size(fVals,1));
    diffs = setdiff(idxs,k(i,:));
    
    %calculate orientation vector
    orientationVector = fVals(diffs(1),:)-refPoint;
    
    orientation = (orientationVector*normal>0);    %check orientation of facet (either 0 or 1)

    normal = normal*(2*orientation-1); %flip normal vector if it faces in the wrong direction
    
    %if the facet is all negative it faces outwards and is irrelevant for
    %the calcuation of the pareto surface
    if all(normal<0)
        %rejected_facets_normal(i,:) = k(i,:);
        continue
    end

    
    %calculate constant for hyperplane equation (c = n*V)
    %n: normal Vector, V: Vector to point on hyperplane
    c = refPoint*normal; %(used for distance calculation in Step 4)

    %%Calculate lower distal point LDP
    facetWeightVectors = penalties(k(i,:),:);   
    %Lower bound approximation for a convex set n(F-P)0.
    %n is the normal of the hyperplane going through point F
    %This leads to a set of equations nP = nF
    rs = sum(facetWeightVectors.*facetPoints,2);
    LDP= linsolve(facetWeightVectors,rs);



    %%Now for LDPs with okay values calculate distance to convexHull
    %conditions that are checked for a facet:
    %LDP should not be below lower limit of bounding box
    %LDP should idealy not be above nadir point component (has to be eased)
    %(LDP should not be "above" PS?)
%     && (((refPoint-LDP')*normal) > 0) 
    if ~(((refPoint-LDP')*normal) > 0)
        'HI'
    end
    if distVal
        upBound = 1;
    else
        upBound = 10;
    end
    if all(LDP>=zeros(size(LDP)))&& all(LDP<=ones(size(LDP))*10) 
        dist = abs(LDP'*normal-c) ;
        facets(i,:) = k(i,:);
        
        normals(i,:) = normal;
        dists(i,1) = dist; 
    end
end

%Now for facets that have not been rejected: Sort by distances

[A,I] = sort(dists,'descend');
maxdist = dists(I(1));
w = zeros(1,size(penalties,2));

found = false;
i = 1;
while ~found && i <= numel(I)
    accuracy = 2;
    idx = I(i);
    norm = normals(idx,:);
    PsPens = penalties(k(idx,:),:);

    %{
 if all(norm>=0)  && ~all(norm == 0) && ~any(ismember(round(penalties,accuracy),round(norm,accuracy),'rows'))
        w = normals(idx,:);
        distfacet = dists(idx,1); 
        found = true;
        zwvec = w;
    %}
    if (~all(normals(idx,:) == 0))&& ~any(ismember(round(penalties,accuracy),round(norm,accuracy),'rows'))
       
        maxminw0 = [5/30,13/30,10/30];
        maxminw0 = maxminw0/sum(maxminw0);
        w = matRad_maxminVector(PsPens,maxminw0);
        w = matRad_maxminVector(PsPens,w);
        zwvec = w;
        w = w/sqrt((sum(w.^2)));
        if all(round(w,accuracy)>=0)  && ~all(round(w,accuracy) == 0) && ~any(ismember(round(penalties,accuracy),round(w,accuracy),'rows'))
            
            distfacet = dists(idx,1); 
            found = true;
        end
    end
    i = i+1;
end
if dists(k(idx,:),:)
    testVar = 1;
end
%%

figure
trisurf(k,fVals(:,1),fVals(:,2),fVals(:,3),'FaceColor','cyan')
hold on
maxDistFacet = fVals(k(idx,:),:)
fill3(maxDistFacet(:,1),maxDistFacet(:,2),maxDistFacet(:,3),'green')

figure
fill3(PsPens(:,1),PsPens(:,2),PsPens(:,3),'red')
hold on 
%%
linest = [zwvec;[0,0,0]]
%%
%plot3(linest(:,1),linest(:,2),linest(:,3))
scatter3(zwvec(1),zwvec(2),zwvec(3))
%%

w
'a'
