function [optTS, U_nom_bank] = pe_param_synth(pe, n_u, n_k, dt)
% Returns CORA createTestSuite options that encode input shape + PE order.
% Note: createTestSuite still draws random values internally for "randn";
% we make parameters deterministic and let CORA do its draw (not identical
% samples, but same excitation shape).

    mode  = lower(char(getf(pe,'mode','randn')));
    L     = max(1, round(getf(pe,'order',2)));
    Aamp  = getf(pe,'strength',1);

    optTS = struct();
    U_nom_bank = []; % reserved for future use 

    switch mode
        case 'randn'
            % CORA's randn curve uses params p=[delay; A; tZero].
            % Choose a window long enough to help achieve PE(L).
            delay = 0;
            tZero = min(n_k, max(ceil(0.3*n_k), ceil(0.5*n_k + L)));
            p = repmat([delay; Aamp; tZero], 1, n_u);

            optTS.inputCurve      = repmat("randn", n_u, 1);
            optTS.inputParameters = p;
            optTS.contInput       = true;

        case 'sinwave'
            % CORA sinWave uses p=[delay; A; T_index; Ts_index], with
            % T = 10*max(dt, p(3))*dt  (per CORA's aux_createCurve)
            % Choose an effective period vs L:
            T_eff = max(8, round(n_k / max(1, L/2))); % in samples
            p3 = T_eff;  p4 = 0;  delay = 0;

            p = repmat([delay; Aamp; p3; p4], 1, n_u);

            optTS.inputCurve      = repmat("sinWave", n_u, 1);
            optTS.inputParameters = p;
            optTS.contInput       = true;

        otherwise
            % Fallback: plain Gaussian
            optTS.inputCurve      = repmat("randn", n_u, 1);
            delay = 0; tZero = min(n_k, max(ceil(0.3*n_k), ceil(0.5*n_k + L)));
            p = repmat([delay; Aamp; tZero], 1, n_u);
            optTS.inputParameters = p;
            optTS.contInput       = true;
    end
end

function v = getf(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
