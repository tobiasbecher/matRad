function orderedPoints = matRad_AdjustedTravellingSalesman(penPoints)
% matRad helper function that solves the Travelling Salesman Problem for a given penalty grid to minimize
% the optimization time for Pareto Point generation. Since no actual loop
% is required, a dummy node is used that connects to all nodes with
% weight 0.
%
% input
%   penPoints:          matrix containing the penalty Points to be reordered
%
% output
%   orderedPoints:      reordered matrix
%
% References
%   https://de.mathworks.com/help/optim/ug/travelling-salesman-problem.html
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
nStops = size(penPoints,1);

%create all possible edges
idxs = nchoosek(1:nStops,2);
nRealEdges = size(idxs,1); % number of edges without dummy node

%create dummy node and create edges to all other nodes
idx2 = zeros(nStops,2);
idx2(1:nStops,1) = nStops + 1;
idx2(1:nStops,2) = 1:nStops;
idxs = [idxs;idx2];

% Calculate distance between all nodes (except for dummy node)
dist = zeros(size(idxs,1),1);
for i = 1:nRealEdges
    dist(i,1) = sqrt(sum((penPoints(idxs(i,1),:)-penPoints(idxs(i,2),:)).^2));
end
lendist = length(dist);

%%Solve Travelling Salesman using integer programming

%create constraints
nStopsMod = nStops + 1;
Aeq = spalloc(nStopsMod,length(idxs),nStopsMod*(nStopsMod-1)); % Allocate a sparse matrix
for ii = 1:nStopsMod
    whichIdxs = (idxs == ii); % Find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % Include trips where ii is at either end
    Aeq(ii,:) = whichIdxs'; % Include in the constint matrix
end
beq = 2*ones(nStopsMod,1);


intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);
%%

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);


x_tsp = logical(x_tsp);
G = graph(idxs(:,1),idxs(:,2),dist);
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2),[],numnodes(G));
%%

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
fprintf('# of subtours: %d\n',numtours);
%%
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
    %hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    %highlight(hGraph,Gsol,'LineStyle','-','EdgeColor','m')
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    fprintf('# of subtours: %d\n',numtours)
end

%use given path to reorder penalty Grid

order = zeros(nStops,1); %initialize ordering
path_edges = [idxs(x_tsp,1),idxs(x_tsp,2)]; % edges along salesman path
% start from edge connected to dummy node 
nodeNum = sum(path_edges(end-1,:))-nStops-1;
order(1,1) = nodeNum;
path_edges(end-1:end,:) = []; % remove dummy node edges

for i = 2:nStops
    [row,col] = find(path_edges == nodeNum);%find connecting edge (a bit ugly)
    nodeNum = sum(path_edges(row,:))-nodeNum;
    order(i,1) = nodeNum;
    path_edges(row,:) = []; % remove dummy node edges
end
orderedPoints = penPoints(order,:);