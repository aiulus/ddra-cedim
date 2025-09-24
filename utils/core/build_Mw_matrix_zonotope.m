function MZ = build_Mw_matrix_zonotope(W, dim_x, T)
%BUILD_MW_MATRIX_ZONOTOPE  Toeplitz lift of a (state-)disturbance zonotope W across T steps.
% Center is replicated; generators are shifted across columns.
   % ---- center: replicate actual W center across T columns ----
   Cw = center(W);                       % [dim_x x 1]
   if isempty(Cw), Cw = zeros(dim_x,1); end
   C  = repmat(Cw, 1, T);                % [dim_x x T]

   % ---- generators: build 3-D array with T shifts per generator ----
   Gw = generators(W);                   % [dim_x x nG]
   nG = size(Gw, 2);
   if nG == 0
       G3 = zeros(dim_x, T, 0);
   else
       G3 = zeros(dim_x, T, nG*T);
       h  = 0;
       for i = 1:nG
           g = Gw(:,i);                  % one generator vector
           firstCol = [g, zeros(dim_x, T-1)];
           for j = 0:T-1
               h = h + 1;
               G3(:,:,h) = circshift(firstCol, [0, j]);
           end
       end
   end

   MZ = matZonotope(C, G3);
end