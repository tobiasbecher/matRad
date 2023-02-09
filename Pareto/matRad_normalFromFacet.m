function [facetPoints,refPoint,normal] = matRad_normalFromFacet(fVals,k,i)
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
    
    %calculate orientation vector
    orientationVector = mean(fVals(unique(k),:))-refPoint;
    
    orientation = (orientationVector*normal>0);    %check orientation of facet (either 0 or 1)

    normal = normal*(2*orientation-1); %flip normal vector if it faces in the wrong direction
    