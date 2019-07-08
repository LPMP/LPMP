# Getting Started

* [Introduction](#introduction)
* [Lagrange decomposition](#lagrange-decomposition)
* [Parser](#parser)
* [Factors](#factors)
* [Message](#message)
* [Problem constructor](#problem-constructor)
* [Solver executable](#solver-executable)
* [Compilation and running](#compilation-and-running)

## Introduction

The LPMP project is a C++ library that facilitates writing efficient Lagrange (dual) decomposition algorithms for integer linear programs (ILP).
Also several solvers for particular ILP classes are provided.
Writing new Lagrange decomposition based solvers in LPMP is often faster than writing them from scratch, since common code that is needed across all solvers is already implemented in such a way that only the problem decomposition specific code needs to be provided additionally. A range of efficient methods for updating Lagrange multipliers is already provided, among them dual block coordinate ascent (a.k.a. message passing) and a proximal bundle method. 
Additionally, new problem decompositions can often utilize existing code stemming from related problems.

Writing a new solver boils down to the following steps:
* Writing a parser for the problem input.
* Writing factors and messages corresponding to the subproblems and Lagrange multipliers of the problem decomposition.
* Writing a ''problem constructor'' taking the problem input and translating it into the problem decomposition. Additionally, the problem constructor usually provides a rounding mechanism to decode primal solutions from any dual solution and possibly cutting plane routines for tightening the LP relaxation provided by the Lagrange decomposition.
* Writing a ''factor message connection''. This connection is a structure linking factors, messages and the problem constructor together. Is is given to the LPMP library which generates a corresponding solver.

We will exemplify writing a simple message passing solver for the quadratic pseudoboolean binary optimization problem (QPBO, a.k.a. QUBO, a.k.a. binary pairwise Markov Random Field) with LPMP.
The QPBO problem can be stated as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{x\in\{0,1\}^n}&space;x^\top&space;Q&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{x\in\{0,1\}^n}&space;x^\top&space;Q&space;x" title="\min_{x\in\{0,1\}^n} x^\top Q x" /></a>.

Note that the solver we write will not be state-of-the-art. There are much more efficient ways to solve LP relaxations equivalent to the one we will optimize, foremost [max-flow based ones](http://pub.ist.ac.at/~vnk/papers/EXTENDED-ROOF-DUALITY.html).

All code lives in the namespace LPMP. The complete example code can be found [here](/doc).

## Lagrange decomposition

We will decompose the QPBO problem into two types of factors corresponding to diagonal entries and off-diagonal ones.
Every diagonal factor corresponds to variable <a href="https://www.codecogs.com/eqnedit.php?latex=x_i&space;=&space;x_i&space;\cdot&space;x_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i&space;=&space;x_i&space;\cdot&space;x_i" title="x_i = x_i \cdot x_i" /></a> and can take two values.
We write the factor as 
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_i&space;\in&space;\{0,1\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_i&space;\in&space;\{0,1\}" title="\mu_i \in \{0,1\}" /></a>.
Every off-diagonal entry corresponds to the product 
<a href="https://www.codecogs.com/eqnedit.php?latex=x_i&space;\cdot&space;x_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i&space;\cdot&space;x_j" title="x_i \cdot x_j" /></a>.
Every off-diagonal factor will hold a tuple of variables 
<a href="https://www.codecogs.com/eqnedit.php?latex=(x_i,x_j,&space;y_{ij})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(x_i,x_j,&space;y_{ij})" title="(x_i,x_j, y_{ij})" /></a>.
We write the corresponding subproblem as
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}&space;=&space;(\mu_{ij}(x_i),&space;\mu_{ij}(x_j),&space;\mu_{ij}(y_{ij}))&space;\in&space;\text{conv}\{&space;(0,0,0),&space;(1,0,0),&space;(0,1,0),&space;(1,1,1)&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}&space;=&space;(\mu_{ij}(x_i),&space;\mu_{ij}(x_j),&space;\mu_{ij}(y_{ij}))&space;\in&space;\text{conv}\{&space;(0,0,0),&space;(1,0,0),&space;(0,1,0),&space;(1,1,1)&space;\}" title="\mu_{ij} = (\mu_{ij}(x_i), \mu_{ij}(x_j), \mu_{ij}(y_{ij})) \in \text{conv}\{ (0,0,0), (1,0,0), (0,1,0), (1,1,1) \}" /></a>.
Hence, in all above combinations for binary <a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}" title="\mu_{ij}" /></a> it holds that <a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}(y_{ij})&space;=&space;\mu_{ij}(x_i)&space;\cdot&space;\mu_{ij}(x_j)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}(y_{ij})&space;=&space;\mu_{ij}(x_i)&space;\cdot&space;\mu_{ij}(x_j)" title="\mu_{ij}(y_{ij}) = \mu_{ij}(x_i) \cdot \mu_{ij}(x_j)" /></a>.
Next, we couple together every off-diagonal factor with the two diagonal factors.
This amounts to
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}(x_i)&space;=&space;\mu_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}(x_i)&space;=&space;\mu_i" title="\mu_{ij}(x_i) = \mu_i" /></a>
and
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}(x_h)&space;=&space;\mu_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}(x_j)&space;=&space;\mu_j" title="\mu_{ij}(x_j) = \mu_j" /></a>

## Parser

First, we define a simple class holding the QPBO problem instance.

```c++
using QPBO_instance = std::vector<std::vector<double>>;
```

Next, given a matrix given in the CSV file format, we write a parser that reads the given file and returns a QPBO_instance.
All the solvers in the LPMP project use the [PEGTL](https://github.com/taocpp/PEGTL) library for parsing files. We do so as well here. We also assume for simplicity that all input numbers are integral.

``` c++
struct number : pegtl::seq< pegtl::opt< pegtl::string<'-'> >, pegtl::star< pegtl::blank >, pegtl::plus< pegtl::digit > > {};

struct entry : pegtl::seq< pegtl::star< pegtl::blank >, number, pegtl::star< pegtl::blank >, pegtl::string<','> > {};
struct row_begin : pegtl::seq<> {};
struct row : pegtl::seq< row_begin, pegtl::star< entry >, pegtl::star< pegtl::blank >, number, pegtl::star< pegtl::blank >, pegtl::eolf > {};

struct grammar : pegtl::plus< row > {};

template< typename Rule >
    struct action
    : pegtl::nothing< Rule > {};

template<> struct action< row_begin > { 
    template<typename Input>
        static void apply(const Input& in, QPBO_instance& input) 
        {   
            input.push_back( {} );
        }   
};  

template<> struct action< number > { 
    template<typename Input>
        static void apply(const Input& in, QPBO_instance& input) 
        {   
            input.back().push_back( std::stod(in.string()) );
        }   
};  

QPBO_instance parse_QPBO_file(const std::string filename)
{
    QPBO_instance q;
    pegtl::file_parser problem(filename);

    const bool read_success = problem.parse<grammar, action>( q );
    return q; 
}
```

A few predefined parser utility classes can be found in [pegtl_parse_rules.h](/include/pegtl_parse_rules.h).

## Factors

Every factor stores the cost of its corresponding subproblem and its primal solution.
It has to implement the following functions: 
* LowerBound(): Given the current cost of the factor, return the value of the optimal solution.
* EvaluatePrimal(): Given a primal solution and costs of the factor, return the cost of the primal assignment,
* init_primal(): Set the values of primal variables to an initial state.
* export_variables(): Expose the cost of the factor to LPMP.
Below, we show the two classes corresponding to diagonal and off-diagonal factors.

```c++
class QPBO_diagonal_factor {
public:
QPBO_diagonal_factor(const double c) : cost(c) {}
double LowerBound() const { return std::min(cost, 0.0); }
double EvaluatePrimal() const { return primal*cost; }
void init_primal() { primal = 0; }
auto export_variables() { return std::tie(cost); }

double cost;
unsigned char primal;
};

class QPBO_off_diagonal_factor {
public:
QPBO_off_diagonal_factor(const double c) : cost({0.0, 0.0, c}) {}
double LowerBound() const {
return std::min({ 0.0, cost[0], cost[1], cost[0] + cost[1] + cost[2] });
}
double EvaluatePrimal() const {
if( primal[0] * primal[1] != primal[2] )
    return std::numeric_limits<double>::infinity();
return primal[0]*cost[0] + primal[1]*cost[1] + primal[2]*cost[2];
}
void init_primal() { primal = {0,0,0}; }
auto export_variables() { return std::tie(cost); }

array<double,3> cost;
std::array<unsigned char,3> primal;
};
```

## Message

A message connects two factors with each other. In LPMP, one factor is called the left one and the other the right one.
In our example, the diagonal factor will be the left, will the off-diagonal factor will be the right one.
Each message implements the following functions: RepamLeft, RepamRight, send_message_to_left, send_message_to_right.
The message will couple diagonal and off-diagonal factors as described in the [Lagrange decomposition](#lagrange-decomposition). We write two classes corresponding to the coupling constraints
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}(x_i)&space;=&space;\mu_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}(x_i)&space;=&space;\mu_i" title="\mu_{ij}(x_i) = \mu_i" /></a>
and
<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{ij}(x_h)&space;=&space;\mu_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{ij}(x_j)&space;=&space;\mu_j" title="\mu_{ij}(x_j) = \mu_j" /></a>
respectively.

```c++
class QPBO_message {
public:
QPBO_message(const unsigned char _index) : index(_index) {}
void RepamLeft(QPBO_diagonal_factor& d, const double val, const std::size_t idx) { d.cost += val; }
void RepamRight(QPBO_off_diagonal_factor& o, const double val, const std::size_t idx) { o.cost[index] += val; }
template<typename MSG>
void send_message_to_right(const QPBO_diagonal_factor& d, MSG& msg, const double scaling)
{
    msg[0] -= scaling*d.cost;
}
    template<typename MSG>
void send_message_to_left(const QPBO_off_diagonal_factor& o, MSG& msg, const double scaling)
{
    const double min_zero_marginal = std::min(0.0, o.cost[1-index]);
    const double min_one_marginal = std::min(o.cost[index], o.cost[0] + o.cost[1] + o.cost[2]);
    msg[0] -= scaling*(min_one_marginal - min_zero_marginal);
}

private:
const unsigned char index;
};

```

For an introduction to the theory behind message passing, in particular what effect the message passing functions have, we refer to this [article](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Dual_Ascent_CVPR_2017_paper.html).

## Problem constructor

The next step is to provide a class taking as input a QPBO_instance and constructing the Lagrange decomposition we outlined above.
The constructor takes a command line object to which it can pass additional parameters and is also given a reference to a solver object. Through the solver object we can factors and messages.
```c++
template<typename FMC, typename DIAGONAL_FACTOR_CONTAINER, typename OFF_DIAGONAL_FACTOR_CONTAINER, typename QPBO_MESSAGE>
class QPBO_problem_constructor {
public:

template<typename SOLVER>
QPBO_problem_constructor(SOLVER& s) : lp(&s.GetLP()) {}

void construct (const QPBO_instance& q)
{
    const std::size_t n = q.size();
    std::vector<DIAGONAL_FACTOR_CONTAINER*> diagonal_factors;
    for(std::size_t i=0; i<n; ++i)
        diagonal_factors.push_back( lp->template add_factor<DIAGONAL_FACTOR_CONTAINER>(q[i][i]) );
    for(std::size_t i=0; i<n; ++i) {
        for(std::size_t j=i+1; j<n; ++j) {
            auto* o = lp->template add_factor<OFF_DIAGONAL_FACTOR_CONTAINER>(q[i][j] + q[j][i]);
            lp->template add_message<QPBO_MESSAGE>(diagonal_factors[i], o, 0);
            lp->template add_message<QPBO_MESSAGE>(diagonal_factors[j], o, 1);
        } 
    }
}

private:
LP<FMC>* lp;
};
```

## Factor message connection

We must finally add a structure that details how factors, messages and the problem constructor fit together. This is done in the ''factor message connection''. To let LPMP know about this relationship, we wrap the factors and messages into a `FactorContainer` and `MessageContainer`' respectively.

```c++
struct QPBO_FMC {
    constexpr static const char* name = "QPBO example"; 

    using diagonal_factor_container = FactorContainer<QPBO_diagonal_factor, QPBO_FMC, 0>;
    using off_diagonal_factor_container = FactorContainer<QPBO_off_diagonal_factor, QPBO_FMC, 1>;

    using QPBO_message_container = MessageContainer<QPBO_message, 0, 1, message_passing_schedule::left, variableMessageNumber, 2, QPBO_FMC, 0>;

    using FactorList = meta::list<diagonal_factor_container, off_diagonal_factor_container>;
    using MessageList = meta::list<QPBO_message_container>;

    using problem_constructor = QPBO_problem_constructor<QPBO_FMC, diagonal_factor_container, off_diagonal_factor_container, QPBO_message_container>;

};
```

For every factor, the template parameters in `FactorContainer` are as follows:
1. The factor class that is used.
2. The ''factor message connection'' which is currently defined.
3. The factor number, consecutively numbered from 0 up.

For every message, the template parameters in `MessageContainer` are as follows:
1. The message class that is used.
2. The factor number of the left factor (the third parameter in the respective  `FactorContainer`).
3. The factor number of the right factor.
4. An identifier that tells LPMP how message updates are scheduled. Possible values are
   * `message_passing_schedule::left`: Whenever left factor is visited, first receive message from right and then send message to right factor.
   * `message_passing_schedule::right`: Whenever right factor is visited, first receive message from left and then send message to left factor.
   * `message_passing_schedule::both`: A combination of `message_passing_schedule::left` and `message_passing_schedule::right`.
   * `message_passing_schedule::none`: Do not update the message.
5. The number of messages that the left factor is connected to of the current type. Possible values are either a positive number, if exactly that number of messages is connected to the respective factor, `variableMessageNumber` if that number is unknown, or `atMostOneMessage`, `atMostTwoMessages`, ... if the number of messages is bounded by the respective number.
6. The number of messages that the right factor is connected to of the current type.
7. The ''factor message message connection'' currently defined.
8. The message number, consecutively numbered from 0 up.

Next, we need to collect all factors in ``FactorList``. Factors must be given in the order of their factor number.
The same mustb e done for messages in ``MessageList``.

Last, the problem constructor is given in ``problem_constructor``.
  

## Solver executable

Finally, we can combine everything together to obtain an executable that solves QPBO instances. To that end we need to provide a visitor. The visitor controls how many iterations are performed, in which intervals primal solutions are decoded, how often cutting plane procedures are invoked etc.

```c++
#include "QPBO_example.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
    Solver<LP<QPBO_FMC>,StandardVisitor> solver(argc,argv);
    auto input = parse_QPBO_file(solver.get_input_file());
    solver.GetProblemConstructor().construct(input);
    return solver.Solve(); 
}
```

## Compilation and running

LPMP uses cmake for the build system. The relevant lines in the CMakeLists.txt file are

```
add_executable(QPBO_example QPBO_example.cpp)
target_link_libraries(QPBO_example LPMP)
```

Most code using LPMP should link against the `LPMP` target so that include paths are set correctly.
Finally, we can run the solver we have written from the command line by typing `./QPBO_example -i ${input_file}`.
