test = matRad_PriorityClass();

test.addObjective(3,DoseObjectives.matRad_SquaredDeviation(1,50),2,1);
test.addObjective(2,DoseObjectives.matRad_SquaredDeviation(1,50),3,1);
test.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),4,1);
test.addConstraint(DoseConstraints.matRad_MinMaxDose(40,50),1)
test.addConstraint(DoseConstraints.matRad_MinMaxDose(40,50),1)
test.addConstraint(DoseConstraints.matRad_MinMaxDose(40,50),2)

load('TG119.mat')
%%
test.generateBasecst(cst)
%%
[cstIt,priority] = test.modifyCst(cst)

%%
test.updateStep(cstIt,'a')


