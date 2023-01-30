function [f,fgrad] =  testDif(al,pen)
            x= al*pen';
            f = - (pen-x).^2';
            fgrad =  (pen'*2*(pen-x))
