results_dir = '';
[T, PS] = load_augmented(results_dir);
audit_augmented(T, PS);
% Core plots
plot_perstep_coverage(PS);
coverage_size_pareto(T);
plot_first_violation(T);
boxplot_directional_error(T);
% Paired tests
paired_compare(T,'E1');
paired_compare(T,'E2');
paired_compare(T,'E3');
% Export compact tables
export_summary_tables(T, fullfile(results_dir,'analysis'));
% Check against artifacts
sanity_check_row(results_dir, 1);
