function [size_scalar, wid_k] = normalize_widths(wid_k, q, agg)
    if nargin<3 || isempty(agg), agg = "mean"; end  
    switch string(agg)
        case "mean"
            wid_k = wid_k ./ q;
        case "sum"
            % no change
        otherwise
            error('normalize_widths: unknown agg=%s', agg);
    end
    size_scalar = mean(wid_k(:), 'omitnan');
end
