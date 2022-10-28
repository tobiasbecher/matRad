[a,b] = matRad_generatePlanarPenaltyGrid(200,[2,2]);
%%
figure
matRad_plotPenaltyGrid(b)
        %%
[a,b] = matRad_generateSphericalPenaltyGrid(200,[1,2]);
%%
figure
matRad_plotPenaltyGrid(b)
%%
c = b./sum(b,2)

%%
matRad_plotPenaltyGrid(c)
%%

