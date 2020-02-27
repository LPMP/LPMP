#pragma once

#include "ILP_input.h"
#include "cuddObj.hh"

namespace LPMP {

    class bdd_solver_interface {

        virtual void init(const ILP_input& instance) = 0;
        virtual double lower_bound() = 0;
        virtual void iteration() = 0;

        // virtual std::size_t nr_bdds() const = 0;
        // virtual std::size_t nr_variables() const = 0;

    };

}
