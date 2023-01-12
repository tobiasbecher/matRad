function totDif = matRad_vectorDiff(al,w)
% matRad_function that calculates the total difference between a vector and
% a matrix of vectors (stored as row vectors)
    wi = (al*w);
    wDif = (w-wi).^2; %subtract vector from each row
    %wDif = sqrt(sum(wDif,2));
    wDif = sum(wDif,2);
    totDif = sum(wDif,1);
end
