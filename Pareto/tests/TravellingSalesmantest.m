%%
nStops = 10;
[a,b] = matRad_generatePlanarPenaltyGrid(nStops,[1,2]);
c = b
%%
idxs = nchoosek(1:nStops,2)
nRealEdges = size(idxs,1)
idx2 = zeros(nStops,2)
idx2(1:nStops,1) = nStops + 1
idx2(1:nStops,2) = 1:nStops
idxs = [idxs;idx2];
%%
dist = zeros(size(idxs,1),1);
for i = 1:nRealEdges
    i
    dist(i,1) = sqrt(sum((b(idxs(i,1),:)-b(idxs(i,2),:)).^2));
end

lendist = length(dist);
%%
'A'
nStopsMod = nStops + 1;
Aeq = spalloc(nStopsMod,length(idxs),nStopsMod*(nStopsMod-1)); % Allocate a sparse matrix
for ii = 1:nStopsMod
    whichIdxs = (idxs == ii) % Find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % Include trips where ii is at either end
    Aeq(ii,:) = whichIdxs'; % Include in the constraint matrix
end
beq = 2*ones(nStopsMod,1);
%%
intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);
%%
opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);
%%
x_tsp = logical(x_tsp)
%%
G = graph(idxs(:,1),idxs(:,2),dist);
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2),[],numnodes(G));

%%
idxs
%%
idxs(x_tsp,2)
%%
Gsol
%%
aaaaaa
stopsLon = b(:,1)
stopsLat = b(:,2)
stopsHor = b(:,3)
stopsLon = [stopsLon;1.0]
stopsLat = [stopsLat;0]
stopsHor = [stopsHor;0]

%%%
figure
hGraph = plot(G,'XData',stopsLon,'YData',stopsLat,'LineStyle','none','NodeLabel',{});
hold on
highlight(hGraph,Gsol,'LineStyle','-')
title('Solution with Subtours')
% Draw the outside border
%%
tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
fprintf('# of subtours: %d\n',numtours);
%%
A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,nStops)]; % A guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1) + 1; % Counter for indexing
        subTourIdx = find(tourIdxs == ii); % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
    x_tsp = logical(round(x_tsp));
    Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2),[],numnodes(G));
    % Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2)); % Also works in most cases
    
    % Visualize result
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-','EdgeColor','m')
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    fprintf('# of subtours: %d\n',numtours)
end
%%
%Gsol.Edges
%%
title('Solution with Subtours Eliminated');
%hold off
%%
%%
x_tsp
idxs