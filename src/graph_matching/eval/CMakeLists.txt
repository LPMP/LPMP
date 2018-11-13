set(EVAL_MODULE_PATH "${PROJECT_SOURCE_DIR}/include/eval")

set(WILLOW_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/WILLOW/instances/WILLOW-ObjectClass_v3/")
configure_file(evaluate_WILLOW.py.in evaluate_WILLOW.py)

set(HOTEL_HOUSE_SYNTHETIC_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/hotel_house_synthetic/")
configure_file(evaluate_hotel_house_synthetic.py.in evaluate_hotel_house_synthetic.py)

set(WORMS_MM_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/worms-mm/")
configure_file(evaluate_worms_mm.py.in evaluate_worms_mm.py)

set(HiPPI_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/HiPPI/")
configure_file(evaluate_HiPPI.py.in evaluate_HiPPI.py)

configure_file(generate_plots_and_table.sh.in generate_plots_and_table.sh)
