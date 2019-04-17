set(EVAL_MODULE_PATH "${PROJECT_SOURCE_DIR}/include/eval")
set(MAX_CUT_EVAL_DIRECTORY "${PROJECT_SOURCE_DIR}/src/max_cut/eval/")
set(MQLIB_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/max_cut/eval/MQLib/")
set(ISING_MODELS_INSTANCE_DIRECTORY "${PROJECT_SOURCE_DIR}/src/max_cut/eval/Ising_models/")

configure_file(download_max_cut_datasets.sh.in download_max_cut_datasets.sh)

configure_file(MQLib_evaluation_data.py.in MQLib_evaluation_data.py)
configure_file(run_MQLib_experiments.py.in run_MQLib_experiments.py)
configure_file(evaluate_MQLib_experiments.py.in evaluate_MQLib_experiments.py)

configure_file(Ising_models_evaluation_data.py.in Ising_models_evaluation_data.py)
configure_file(run_Ising_models_experiments.py.in run_Ising_models_experiments.py)
configure_file(evaluate_Ising_models_experiments.py.in evaluate_Ising_models_experiments.py) 
