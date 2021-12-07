% Esteban Vazquez-Hidalgo
% last update 07.13.2021
% mappingMatrix.m creates a matrix of neighboring cells for current
% filament. 
% find number of elements
tots = rows*cols;
% a variable holding the number of rows and columns
ss = [rows, cols];
% creat a vector to hold all the neighbors, where the row number is the
% current filament and the columns of that row are the neighbors of the
% column number
mapps = zeros(tots,8);
% vectur used to find neighbors
Ni = [-1 -1 -1  0  0  1 1 1];
Nj = [-1  0  1 -1  1 -1 0 1];
% gets the index from the linear value
% the input is the linear index, the output is the subcript index
for kk = 1:tots
    [I, J] = ind2sub(ss,kk);
    % this gets the neighbors for I,J
    tempi = I + Ni;
    tempj = J + Nj;
    % this ensures the periodic boundary conditions are met and that no
    % index is out of the range
    tempi(tempi == 0) = rows;
    tempi(tempi > rows) = 1;
    tempj(tempj == 0) = cols;
    tempj(tempj > cols) = 1;
    % create a temporary vector for each row
    for k = 1:8
        tempr(k) = [tempi(k)];
        tempc(k) = [tempj(k)];
    end
    % this creates a temporary vector that hold the subscripts of the
    % neighbors
    tempv = [tempr' tempc'];
    % this converts the subscripts back into linear indeces. use these linear
    % indices to check neighboring attachment status
    for k = 1:8
        ind(k) = sub2ind(ss, tempv(k,1),tempv(k,2));
    end
    bbb = ind;
    mapps(kk,:) = ind;
end



