#include "lifted_disjoint_paths/lifted_disjoint_paths_input.h"



// implementation of file reading

namespace LPMP {

    namespace lifted_disjoint_paths {


    LdpInstance<> parse_file(const std::string& filename){
    	ConfigDisjoint<> configParams(filename);
    	LdpInstance<> ldpInstance(configParams);
    	return ldpInstance;
    }



    }
}
