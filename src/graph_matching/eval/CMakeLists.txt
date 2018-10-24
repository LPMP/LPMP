set(WILLOW_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/WILLOW/instances/WILLOW-ObjectClass_v3/")
configure_file(evaluate_WILLOW.py.in evaluate_WILLOW.py)

set(HOTEL_HOUSE_SYNTHETIC_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/graph_matching/eval/hotel_house_synthetic/")
configure_file(evaluate_hotel_house_synthetic.py.in evaluate_hotel_house_synthetic.py)
