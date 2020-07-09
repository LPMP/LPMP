set(EVAL_MODULE_PATH "${PROJECT_SOURCE_DIR}/src/bdd/eval")

set(GRAPH_MATCHING_INSTANCE_DIRECTORY "/BS/discrete_opt/work/datasets/graph_matching/")
configure_file(data_graph_matching.py.in data_graph_matching.py)
configure_file(run_graph_matching_experiments.py.in run_graph_matching_experiments.py)
configure_file(evaluate_graph_matching_experiments.py.in evaluate_graph_matching_experiments.py)

set(CELL_TRACKING_INSTANCE_DIRECTORY "/BS/discrete_opt/work/datasets/cell-tracking/")
configure_file(data_cell_tracking.py.in data_cell_tracking.py)
configure_file(run_cell_tracking_experiments.py.in run_cell_tracking_experiments.py)
configure_file(evaluate_cell_tracking_experiments.py.in evaluate_cell_tracking_experiments.py)

set(MRF_INSTANCE_DIRECTORY "/BS/discrete_opt/work/datasets/mrf/")
configure_file(data_mrf.py.in data_mrf.py)
configure_file(run_mrf_experiments.py.in run_mrf_experiments.py)
configure_file(evaluate_mrf_experiments.py.in evaluate_mrf_experiments.py)

set(DISCRETE_TOMOGRAPHY_INSTANCE_DIRECTORY "/BS/discrete_opt/work/datasets/discrete_tomography/discrete_tomography_synthetic/mp/")
configure_file(data_discrete_tomography.py.in data_discrete_tomography.py)
configure_file(run_discrete_tomography_experiments.py.in run_discrete_tomography_experiments.py)
configure_file(evaluate_discrete_tomography_experiments.py.in evaluate_discrete_tomography_experiments.py)

#configure_file(generate_plots_and_table.sh.in generate_plots_and_table.sh)
