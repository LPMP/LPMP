#ifndef LPMP_HORIZON_TRACKING_UAI_INPUT_H
#define LPMP_HORIZON_TRACKING_UAI_INPUT_H

#include <vector>
#include "mrf/mrf_input.h"

namespace LPMP {

    struct horizon_tracking_input {
        mrf_input mrf;
        std::vector<mrf_input> bottleneck_potentials;
    };

    namespace horizon_tracking_uai_input {

        horizon_tracking_input parse_file(const std::string& filename);
        horizon_tracking_input parse_string(const std::string& filename);

    } // namespace horizon_tracking_input

} // namespace LPMP 

#endif // LPMP_HORIZON_TRACKING_UAI_INPUT_H
