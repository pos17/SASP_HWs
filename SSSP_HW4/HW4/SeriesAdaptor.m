function portMatrix =  SeriesAdaptor(Z1,Z2,Z3,adaptCoeff)
    portMatrix = zeros(3);
    totZ = sum([Z1 Z2 Z3]);
    onesRow = ones(1,3);
    onesDiag = diag(onesRow);
