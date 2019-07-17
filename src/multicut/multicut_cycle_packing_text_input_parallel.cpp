#include "multicut/multicut_cycle_packing_parallel.h"
#include "multicut/multicut_text_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
    if(argc != 3)
        throw std::runtime_error("[prog_name] [input_file] [nr_threads]");
    auto input = LPMP::multicut_text_input::parse_file(argv[1]);
    const int nr_thread = std::stoi(argv[2]);
    multicut_cycle_packing_parallel(input, nr_thread);
}
