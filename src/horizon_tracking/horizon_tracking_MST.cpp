#include "mrf/mrf_uai_input.h"
#include "horizon_tracking/horizon_tracking_MST_solver.hxx"

using namespace LPMP;

int main(int argc, char** argv) {
    TCLAP::CmdLine cmd("Horizon Tracking using Minimum Spanning Tree");
    TCLAP::ValueArg<std::string> inputFileArg("i","inputFile","file from which to read problem instance",false,"","file name",cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o","outputFile","file to write solution",false,"","file name",cmd);
    cmd.parse(argc,argv);
    std::string inputFile(inputFileArg.getValue());
    std::string outputFile(outputFileArg.getValue());

    auto input = mrf_uai_input::parse_file(inputFile);
    auto mst_solver = horizon_tracking_MST_solver(input);
    mst_solver.ComputeSolution();
    mst_solver.PrintPrimal();
    if(outputFileArg.isSet()) {
        mst_solver.WritePrimal(outputFile);
    }
}