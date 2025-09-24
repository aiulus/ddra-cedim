function Z2 = linearMap(Z, A)
% compatibility shim when CORA exposes affineMap instead of linearMap
if exist('affineMap','file')
    Z2 = affineMap(Z, A);
else
    error('linearMap/affineMap not found on path. Check CORA installation.');
end
end

