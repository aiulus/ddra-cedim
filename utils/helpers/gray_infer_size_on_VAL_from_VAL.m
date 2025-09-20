function [sizeI_gray, wid_gray_k] = gray_infer_size_on_VAL_from_VAL(sys, VAL, C, params, varargin)
    nk = C.shared.n_k_val;
    wid_all = zeros(nk, numel(VAL.y));
    for i = 1:numel(VAL.y)
        tc = testCase(reshape(VAL.y{i}, nk, [], 1), VAL.u{i}, VAL.x0{i}, sys.dt, 'linearSysDT');
        [~, wk] = gray_infer_size_on_VAL(sys, {tc}, C, params, varargin{:});
        wid_all(:,i) = wk(:);
    end
    wid_gray_k = mean(wid_all, 2, 'omitnan');
    sizeI_gray = mean(wid_gray_k, 'omitnan');
end
