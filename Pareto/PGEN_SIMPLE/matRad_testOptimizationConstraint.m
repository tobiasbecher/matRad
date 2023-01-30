function [c,ceq] = matRad_testOptimizationConstraint(x)
    c = [(x(2)-9)^2+(x(3)-3)^2-x(1),(x(1)-4)^2+(x(3)-3)^2-x(2),(x(1)-4)^2+(x(2)-9)^2-x(3)];
    ceq = [];