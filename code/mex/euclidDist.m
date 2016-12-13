function distance = euclidDist(x1,x2)
% calculates the Euclidean distance from x1[1x2] to x2[1x2]
    distance = sqrt(sum((x1-x2).^2));
end %of function euclidDist


