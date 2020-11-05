#ifndef LDP_PATH_SEPARATOR_HXX
#define LDP_PATH_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP {

template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR>
class ldp_path_separator {

public:
    ldp_path_separator(std::vector<PATH_FACTOR*>& _factorContainer, std::vector<PATH_SNC_MESSAGE*>&  _messageContainer, std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& _sncFactorContainer, const lifted_disjoint_paths::LdpInstance * _pInstance):
    factorContainer(_factorContainer),
    messageContainer(_messageContainer),
    pInstance(_pInstance),
    sncFactorContainer(_sncFactorContainer)
    {

    }




private:
std::vector<PATH_FACTOR*>& factorContainer;
std::vector<PATH_SNC_MESSAGE*>&  messageContainer;
const lifted_disjoint_paths::LdpInstance * pInstance;
std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& sncFactorContainer;


};





}



#endif // LDP_PATH_SEPARATOR_HXX
