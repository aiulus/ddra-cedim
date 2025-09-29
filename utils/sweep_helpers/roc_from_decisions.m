function ROC = roc_from_decisions(safe_true, safe_alg_tau)
% safe_true:  (N x 1) logical
% safe_alg_tau.ddra / .gray: (N x T) logical (safe decision per tau)

algos = {'ddra','gray'};
for a = 1:numel(algos)
    A = algos{a};
    P = safe_alg_tau.(A);   % N x T
    N = numel(safe_true);
    T = size(P,2);

    TPR = zeros(T,1); FPR = zeros(T,1);
    TP = zeros(T,1);  FP = zeros(T,1);
    TN = zeros(T,1);  FN = zeros(T,1);

    pos = safe_true;           nPos = max(1,sum(pos));
    neg = ~safe_true;          nNeg = max(1,sum(neg));

    for t = 1:T
        pred = P(:,t);         % predicted "safe"
        TP(t) = sum(pred &  pos);
        FP(t) = sum(pred & ~pos);
        TN(t) = sum(~pred & ~pos);
        FN(t) = sum(~pred &  pos);

        TPR(t) = TP(t)/nPos;
        FPR(t) = FP(t)/nNeg;
    end

    % Sort by FPR to compute AUC
    [FPRs, idx] = sort(FPR);
    TPRs = TPR(idx);
    AUC  = trapz(FPRs, TPRs);

    ROC.(A) = struct('TPR',TPR,'FPR',FPR,'TP',TP,'FP',FP,'TN',TN,'FN',FN,'AUC',AUC);
end
end
