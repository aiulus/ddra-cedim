function h = val_hash(v)
%VAL_HASH Deterministic hash of structs/cells/arrays for parity checks.
% Returns lowercase hex SHA-256 string.

    md = java.security.MessageDigest.getInstance('SHA-256');
    sweep_feed(md, v);
    h = lower(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
end