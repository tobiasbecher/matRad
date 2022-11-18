function [normC,ceq] =  matRad_normConstraint(al,w)
    wi = al*w;
    normC = wi*wi'-1;
    ceq =[];