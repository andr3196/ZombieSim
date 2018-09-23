function  MAsMatrix = toColumnInverseN(M,n)
global dx dy
MAsMatrix = reshape(M,dy, dx, n);