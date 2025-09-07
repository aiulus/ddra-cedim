function env_report
fprintf('--- MATLAB ---\n');
fprintf('MATLAB: %s (%s)\n', version, version('-release'));

opt = ver('optim');
if ~isempty(opt)
    fprintf('Optimization Toolbox: %s (v%s)\n', opt.Name, opt.Version);
else
    fprintf('Optimization Toolbox: NOT FOUND\n');
end

% Try CORA via ver
v = ver;
vCora = v(contains({v.Name},'CORA','IgnoreCase',true) | contains({v.Name},'Reach','IgnoreCase',true));
if ~isempty(vCora)
    fprintf('CORA: %s (v%s)\n', vCora(1).Name, vCora(1).Version);
else
    % Heuristic path discovery
    p = which('zonotope');
    if ~isempty(p)
        coraRoot = fileparts(fileparts(p));
        fprintf('CORA root: %s\n', coraRoot);
        try
            if exist(fullfile(coraRoot,'version.txt'),'file')
                vtxt = strtrim(fileread(fullfile(coraRoot,'version.txt')));
                fprintf('CORA version.txt: %s\n', vtxt);
            end
        catch, end
        % Git info
        [~,tag] = system(sprintf('git -C "%s" describe --tags --dirty --always', coraRoot));
        if ~isempty(strtrim(tag))
            fprintf('CORA git: %s\n', strtrim(tag));
        end
    else
        fprintf('CORA: NOT FOUND ON PATH\n');
    end
end

fprintf('\n--- Platform ---\n');
[arch,~,~] = computer;
fprintf('computer: %s\n', arch);
fprintf('%s\n', system_dependent('getos'));

if ispc
    [~,cpu] = system('wmic cpu get Name /value');
elseif ismac
    [~,cpu] = system('sysctl -n machdep.cpu.brand_string');
else
    [~,cpu] = system('lscpu | grep -m1 "Model name" | cut -d: -f2-');
end
fprintf('CPU: %s\n', strtrim(cpu));
fprintf('Cores: %d\n', feature('numcores'));

if ispc
    try; m = memory; fprintf('RAM (Total/Avail): %.1f/%.1f GB\n', m.MemTotalPhys/2^30, m.MemAvailableAllArrays/2^30); end
elseif ismac
    system('vm_stat | head');
else
    system('free -h');
end
end
