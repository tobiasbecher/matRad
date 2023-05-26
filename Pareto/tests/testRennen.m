[ps,errors] = matRad_Rennen(49)
%%


[k,facets] = matRad_ParetoSurfFromFacets(ps)
%%
'a'
figure
trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor',[0.9,0.9,0.9])
hold on 
scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
        'MarkerFaceColor',[1 1 1])
xlabel('f_1(x)')
ylabel('f_2(x)')
zlabel('f_3(x)')
zlim([0.5,5])
ylim([7,9])
xlim([1.75,5])
5%%
figure
plot(errors)