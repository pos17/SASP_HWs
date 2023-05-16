function portMatrix =  SeriesAdaptor(Z1,Z2,Z3)
    portMatrix = zeros(3);
    totZ = sum([Z1 Z2 Z3]);
    onesRow = ones(1,3);
    onesDiag = diag(onesRow);
    colZ = [Z1;Z2;Z3];
    portMatrix = onesDiag - ((2/totZ)*(colZ)*onesRow);
