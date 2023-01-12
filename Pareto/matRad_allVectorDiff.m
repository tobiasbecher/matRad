function wDif = matRad_allVectorDiff(al,w)
    % matRad_function that calculates the difference between a 
        wi = (al*w)
        wDif = (w-wi).^2 %subtract vector from each row
        wDif = sqrt((sum(wDif,2)))
        wDif = -1*wDif;
    end
    