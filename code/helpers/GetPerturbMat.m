function perturbMat = GetPerturbMat(matType,rad,neigh)
% perturbMat = GetPerturbMat(matType,rad)

if ~exist('rad','var')
    rad = 3;
end
if ~exist('neigh','var')
    neigh = 4;
end
assert(rad==3 || rad==5)
assert(neigh==4 || neigh==8)

if neigh == 8
    perturbMat = zeros(8,rad^8,matType); % 8 x 6561
    
    if rad==3
        vec = -1:1;
    elseif rad==5
        vec = -2:2;
    end
    
    ind = 0;
    for i1 = vec
        for i2 = vec
            for i3 = vec
                for i4 = vec
                    for i5 = vec
                        for i6 = vec
                            for i7 = vec
                                for i8 = vec
                                    ind = ind + 1;
                                    perturbMat(:,ind) = [i1;i2;i3;i4;i5;i6;i7;i8];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if neigh == 4
    perturbMat = zeros(8,5^4,matType); % 8 x 6561
    
    if rad==3
        vec = -1:1;
    elseif rad==5
        error('dont use this config'); % vec = -2:2;
    end
    
    ind = 0;
    for i1 = vec
        for i2 = vec
            if (abs(i1) + abs(i2) == 2)
                continue
            end
            for i3 = vec
                for i4 = vec
                    if (abs(i3) + abs(i4) == 2)
                        continue
                    end
                    for i5 = vec
                        for i6 = vec
                            if (abs(i5) + abs(i6) == 2)
                                continue
                            end
                            for i7 = vec
                                for i8 = vec
                                    if (abs(i7) + abs(i8) == 2)
                                        continue
                                    end
                                    ind = ind + 1;
                                    perturbMat(:,ind) = [i1;i2;i3;i4;i5;i6;i7;i8];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
'';
