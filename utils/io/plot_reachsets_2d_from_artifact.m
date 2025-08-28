function plot_reachsets_2d_from_artifact(artifact_path, b, dims, varargin)
    S = load(artifact_path);  
    plot_reachsets_2d(S.sys_gray, S.sys_ddra, S.VAL, S.W_eff, b, dims, varargin{:});
end
