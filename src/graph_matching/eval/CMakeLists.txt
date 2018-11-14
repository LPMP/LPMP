set(EVAL_MODULE_PATH "${PROJECT_SOURCE_DIR}/include/eval")

set(WILLOW_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/WILLOW/instances/WILLOW-ObjectClass_v3/")
configure_file(WILLOW_evaluation_data.py.in WILLOW_evaluation_data.py)
configure_file(run_WILLOW_experiments.py.in run_WILLOW_experiments.py)
configure_file(evaluate_WILLOW_experiments.py.in evaluate_WILLOW_experiments.py)

set(HOTEL_HOUSE_SYNTHETIC_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/hotel_house_synthetic/")
configure_file(hotel_house_synthetic_evaluation_data.py.in hotel_house_synthetic_evaluation_data.py)
configure_file(run_hotel_house_synthetic_experiments.py.in run_hotel_house_synthetic_experiments.py)
configure_file(evaluate_hotel_house_synthetic_experiments.py.in evaluate_hotel_house_synthetic_experiments.py)

set(WORMS_MM_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/worms-mm/")
configure_file(worms_mm_evaluation_data.py.in worms_mm_evaluation_data.py)
configure_file(run_worms_mm_experiments.py.in run_worms_mm_experiments.py)
configure_file(evaluate_worms_mm_experiments.py.in evaluate_worms_mm_experiments.py)

set(HiPPI_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/HiPPI/")
configure_file(evaluate_HiPPI.py.in evaluate_HiPPI.py)

configure_file(generate_plots_and_table.sh.in generate_plots_and_table.sh)
