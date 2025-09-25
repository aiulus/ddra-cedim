function s = stable_seed(D, alpha_w, pe, dyn, typ, tag)
    o     = double(getfielddef(pe,'order',0));
    smode = double(sum(char(string(getfielddef(pe,'mode','')))));
    dyn   = char(string(dyn));
    typ   = char(string(typ));
    tag   = char(string(tag));
    s64   = uint64(1469598103934665603);
    K     = uint64([D, round(alpha_w*1e6), o, smode, sum(double(dyn)), sum(double(typ)), sum(double(tag))]);
    for v = K
        s64 = bitxor(s64, uint64(v));
        s64 = s64 * 1099511628211;
    end
    s = uint32(mod(s64, uint64(2^31-1))); if s==0, s=1; end
end
