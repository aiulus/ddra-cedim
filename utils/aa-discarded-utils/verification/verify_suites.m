function verify_suites(TS, U_nom_all, n_k)
    assert(iscell(TS));
    for m = 1:numel(TS)
        % u should be 2-D and equal to nominal
        assert(ismatrix(TS{m}.u) && size(TS{m}.u,1)==n_k, 'TS{%d}.u shape wrong', m);
        if norm(TS{m}.u - U_nom_all{m}.', Inf) > 1e-12
            error('TS{%d}.u != U_nom_all{%d}', m, m);
        end
        % y is 3-D with samples in 3rd dim
        assert(ndims(TS{m}.y)==3 && size(TS{m}.y,1)==n_k, 'TS{%d}.y shape wrong', m);
    end
    disp('OK: testSuite structure canonical for gray.');
end
