[ps,errors] = matRad_Rennen(50)
%%


[k,facets] = matRad_ParetoSurfFromFacets(ps)
%%
figure
trisurf(k,ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')

figure
trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')
hold on 
scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
        'MarkerFaceColor',[0 0 0])

%%
figure
plot(errors)