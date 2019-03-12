#include <random>
#include "test.h"
#include "../test_message.hxx"
#include "max_cut/max_cut_factors_messages.h"

using namespace LPMP;
int main(int argc, char** argv)
{
    std::random_device rd;

    max_cut_triplet_factor t012;
    max_cut_triplet_factor t013;
    max_cut_triplet_factor t014;
    max_cut_triplet_factor t023;
    max_cut_triplet_factor t024;
    max_cut_triplet_factor t034;
    max_cut_triplet_factor t123;
    max_cut_triplet_factor t124;
    max_cut_triplet_factor t134;
    max_cut_triplet_factor t234;
    max_cut_quintuplet_factor q;

    max_cut_triplet_quintuplet_message_012::msg_val_type msg_val{};

    max_cut_triplet_quintuplet_message_012 msg_012;
    max_cut_triplet_quintuplet_message_013 msg_013;
    max_cut_triplet_quintuplet_message_014 msg_014;
    max_cut_triplet_quintuplet_message_023 msg_023;
    max_cut_triplet_quintuplet_message_024 msg_024;
    max_cut_triplet_quintuplet_message_034 msg_034;
    max_cut_triplet_quintuplet_message_123 msg_123;
    max_cut_triplet_quintuplet_message_124 msg_124;
    max_cut_triplet_quintuplet_message_134 msg_134;
    max_cut_triplet_quintuplet_message_234 msg_234;

    // TODO: currently there is no MaximizePotentialAndCOmputePrimal method that would work with partial labelings, hence test_repam_message etc. does not work.
    /*
    test_repam_message<Chirality::left>(t012, q, msg_012, msg_val, rd);
    test_repam_message<Chirality::left>(t013, q, msg_013, msg_val, rd);
    test_repam_message<Chirality::left>(t012, q, msg_014, msg_val, rd);
    test_repam_message<Chirality::left>(t023, q, msg_023, msg_val, rd);
    test_repam_message<Chirality::left>(t024, q, msg_024, msg_val, rd);
    test_repam_message<Chirality::left>(t034, q, msg_034, msg_val, rd);
    test_repam_message<Chirality::left>(t123, q, msg_123, msg_val, rd);
    test_repam_message<Chirality::left>(t124, q, msg_124, msg_val, rd);
    test_repam_message<Chirality::left>(t134, q, msg_134, msg_val, rd);
    test_repam_message<Chirality::left>(t234, q, msg_234, msg_val, rd);
    */
}
