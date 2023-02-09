function [ps,errors] = matRad_Rennen(nIter)
    
    %%Initialize OPS matrices
    OPSA = [];
    OPSb = [];


    %% Part should be removed later on    
    %calculate initial points
    penalties = [[1,0,0];[0,1,0];[0,0,1];];
    initSize = size(penalties,1);
    ps = zeros(size(penalties)+[nIter,0]);

    for i = 1:size(penalties,1)
        ps(i,:) = matRad_testOptimization(penalties(i,:),[[0,0,0];[1,1,1]],[0,1,0]);
    end

    %calculate first non anchor point
    np = size(penalties,1);
    [a,b,firstNormal] = matRad_normalFromFacet(ps,1:np,1); %calculate normal vector of initial facet
    penalties(np+1,:) = abs(firstNormal);
    
    ps(np+1,:) = matRad_testOptimization(penalties(np+1,:),[[0,0,0];[1,1,1]],[0,1,0]);
    %
    OPSA = [OPSA;-penalties(np+1,:)];
    OPSb = [OPSb;-1*ps(np+1,:)*penalties(np+1,:)'];


    %% Remaining steps
    penalties = [penalties;zeros(nIter,size(penalties,2))]
    errors = [];


    for i = 1:nIter
        %Step 1 calculate convex Hull -> Inner approximation (IPS) and gives facets
        %Rennen Algorithm
        fVals = ps(1:size(penalties,1)-nIter+i-1,:);
        
        %calculate epsilon value
        L = min(fVals,[],1);
        U = max(fVals,[],1);
        eps = U - L;


        fValsMod = matRad_generateDummyPoints(fVals); %generate dummy points
        %
        [k,vol] = convhulln(fValsMod);
        [kred,vol] = convhulln(fVals);
        %check for relevant facets (those that contain points of the original
        %fVals set)
        IPSidxs = 1:size(fVals,1);
        relFacetidxs = [];
                
        for j = 1:size(k,1)
            if any(ismember(k(j,:),IPSidxs))
                relFacetidxs = [relFacetidxs,j];
            end
        end
        facetMods = k(relFacetidxs,:);
        facetErrors = zeros(size(facetMods,1),1);
        normals = zeros(size(facetMods));
        
        %% plots
        
        %hold on
        %{
        figure
        trisurf(kred,ps(1:i+4-1,1),ps(1:i+4-1,2),ps(1:i+4-1,3),'FaceColor','red') 
        %}
    
        %Loop over facets ands 
        for j = 1:size(relFacetidxs,2)
            [facetPoints,refPoint,normal] = matRad_normalFromFacet(fValsMod,facetMods,j);

            
            %check for sign of normals (should be redundant)
            if all(normal<0)
                continue
            end
            
            %now check for OPS point for facet
            lb = min(fVals,[],1);
            ub = max(fVals,[],1);
            z = linprog(normal,OPSA,OPSb,[],[],lb,ub); 
            
            %hyperplane distance
            b = refPoint*normal;

            %calculate error for each facet
            
            facetErrors(j) = (b-z'*normal)/(eps*normal); %calculation is hopefully correct
            normals(j,:) = normal;
            %{
            figure
            trisurf(k,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
            hold on 
            fill3(facetPoints(:,1),facetPoints(:,2),facetPoints(:,3),'green')
            %}
        end

        [A,I] = sort(facetErrors,'descend');

        %%check for next facet to run
        found = false;
        facetNum= 1;    
        w = zeros(1,size(penalties,2));
        accuracy = 3;

        while ~found && facetNum <= numel(I) %loop over facets
            idx = I(facetNum);
            norm = normals(idx,:);
            if ~any(ismember(round(penalties,accuracy),round(norm,accuracy),'rows'))
                errors = [errors,facetErrors(idx)];
                w = norm;
                %run point 
                penalties(i+4,:) = w;
                ps(i+4,:) = matRad_testOptimization(w,[[0,0,0];[1,1,1]],[0,1,0]); %!hardcoded right now
                found = true;
                %%              
                %{
                figure
                trisurf(k,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
                hold on
                fill3(fValsMod(facetMods(idx,:),1),fValsMod(facetMods(idx,:),2),fValsMod(facetMods(idx,:),3),'r')
                %}
                %%
            end
            facetNum = facetNum +1;
        end    

        
        %run for facet and check whether point has already been calculated or not 
        

        % when final point is found: Update OPsw and OPSb
        OPSA = [OPSA;-w]; %add normal vector of facet that was run 
        OPSb = [OPSb;-ps(i+4,:)*w'];

    end
    
