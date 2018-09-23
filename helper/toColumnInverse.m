function  MAsMatrix = toColumnInverse(M)
global dx dy
MAsMatrix = reshape(M,dy, dx, 4);