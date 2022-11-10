function totDif = matRad_vectorDiff(wi,w)
% matRad_function that calculates the total difference between a vector and
% a matrix of vectors (stored as row vectors)
    wDif = (w-wi).^2; %subtract vector from each row
    wDif
    totDif = sum(wDif,'all');
    

end
