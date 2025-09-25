function IO = sweepio_begin(cfg, plots_dir, results_dir, Ntot)
    if nargin < 4, Ntot = []; end
    IO = struct();
    IO.plots_dir    = plots_dir;
    IO.results_dir  = results_dir;
    IO.csv_path     = fullfile(results_dir, 'summary.csv');
    IO.csv_perstep  = fullfile(results_dir, 'summary_perstep.csv');
    IO.append       = getfielddef(getfielddef(cfg,'lowmem',struct()), 'append_csv', false);
    IO.Ntot         = Ntot;
    IO.rowi         = 0;
    IO.hdr          = {};
    IO.cells        = [];                 
    IO.schema_ver   = 2;                

    if IO.append && exist(IO.csv_path,'file'), delete(IO.csv_path); end
    if exist(IO.csv_perstep,'file'), delete(IO.csv_perstep); end

    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    if ~exist(plots_dir,'dir'),  mkdir(plots_dir);  end
end
