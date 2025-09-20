function VAL_report(A,B,tag)
%VAL_REPORT Quick diagnostic between two VAL structs.
    if nargin<3, tag=''; end
    ha = val_hash(A); hb = val_hash(B);
    fprintf('[VAL_report%s] hash(A)=%s  hash(B)=%s  equal=%d\n', ...
        iff(~isempty(tag), [':' tag], ''), ha, hb, strcmp(ha,hb));

    la = numel(A.y); lb = numel(B.y);
    fprintf('  lens: A=%d  B=%d\n', la, lb);
    if la==0 || lb==0, return; end

    nkA = size(A.y{1},1); qA = size(A.y{1},2);
    nkB = size(B.y{1},1); qB = size(B.y{1},2);
    pA = size(A.u{1},2);  pB = size(B.u{1},2);
    nxA= size(A.x0{1},1); nxB= size(B.x0{1},1);
    fprintf('  dims: nk(A,B)=(%d,%d)  q(A,B)=(%d,%d)  p(A,B)=(%d,%d)  nx(A,B)=(%d,%d)\n', ...
        nkA,nkB,qA,qB,pA,pB,nxA,nxB);

    n = min(la,lb);
    for i=1:n
        d0 = max(abs(A.x0{i}-B.x0{i}),[],'all');
        du = max(abs(A.u{i} -B.u{i} ),[],'all');
        dy = max(abs(A.y{i} -B.y{i} ),[],'all');
        if any([d0,du,dy] > 0)
            fprintf('  first mismatch at i=%d: d0=%.3g du=%.3g dy=%.3g\n', i, d0, du, dy);
            break
        end
    end
end

function y = iff(c,a,b)
    if c, y=a; else, y=b; end 
end 

% -------------------- HELPERS --------------------

% tests/helpers/must.m
function must(cond,msg), if ~cond, error("TEST:fail -> " + string(msg)); end, end

% tests/helpers/hash_bytes.m  (stable SHA-1 over numeric/struct/zono)
function h = hash_bytes(x)
    % Normalizes common types to a byte stream, then SHA-1.
    bytes = to_bytes(x);

    md = java.security.MessageDigest.getInstance('SHA-1');
    md.update(uint8(bytes));
    raw = typecast(md.digest(),'uint8');
    hex = lower(reshape(dec2hex(raw).',1,[]));
    h = string(hex);
end

function bytes = to_bytes(x)
    if isa(x,'zonotope')
        c = center(x);            % <-- temp var
        G = generators(x);        % <-- temp var
        s = struct('t','zono','c',double(c),'G',double(G));
        bytes = getByteStreamFromArray(orderfields(s));

    elseif isnumeric(x) || islogical(x)
        bytes = getByteStreamFromArray(double(x));

    elseif isstring(x)
        bytes = getByteStreamFromArray(char(x));

    elseif ischar(x)
        bytes = getByteStreamFromArray(x);

    elseif isstruct(x)
        bytes = getByteStreamFromArray(orderfields(x));

    else
        % Fallback: best-effort serialization
        bytes = getByteStreamFromArray(x);
    end
end


% tests/helpers/val_hash.m (VAL struct -> single hash)
function h = val_hash(VAL)
    parts = strings(0,1);
    for b = 1:numel(VAL.x0)
        parts(end+1) = hash_bytes(VAL.x0{b});
        parts(end+1) = hash_bytes(VAL.u{b});
        parts(end+1) = hash_bytes(VAL.y{b});
    end
    h = hash_bytes(char(join(parts,"|")));
end

% tests/helpers/interval_contains_interval.m
function tf = interval_contains_interval(A, B, tol)
    if nargin<3, tol=1e-9; end
    IA = interval(A); IB = interval(B);
    try, alo=infimum(IA); ahi=supremum(IA); catch, alo=IA.inf; ahi=IA.sup; end
    try, blo=infimum(IB); bhi=supremum(IB); catch, blo=IB.inf; bhi=IB.sup; end
    tf = all(blo >= alo - tol) && all(bhi <= ahi + tol);
end
