function K = red_order(C)
% Unified reduction order:
%   primary  : C.shared.options_reach.zonotopeOrder
%   fallback : 100 (keeps old behavior sane if missing)
    K = 100;
    if isfield(C,'shared') && isfield(C.shared,'options_reach') ...
            && isfield(C.shared.options_reach,'zonotopeOrder') ...
            && ~isempty(C.shared.options_reach.zonotopeOrder)
        K = C.shared.options_reach.zonotopeOrder;
    end
    % Optional compatibility cap:
    if isfield(C,'lowmem') && isfield(C.lowmem,'zonotopeOrder_cap') ...
            && ~isempty(C.lowmem.zonotopeOrder_cap)
        K = min(K, C.lowmem.zonotopeOrder_cap);
    end
end
