#ifndef LPMP_STANDARD_VISITOR_HXX
#define LPMP_STANDARD_VISITOR_HXX

#include "LP.h"
#include "config.hxx"
#include "mem_use.c"
#include "tclap/CmdLine.h"
#include <chrono>

/*
 minimal visitor class:
class Visitor {
public:
   Visitor(TCLAP::CmdLine& cmd);
   LpControl begin(LP& lp);
   LpControl visit(LpControl, const REAL lower_bound, const REAL primal)
};
*/

namespace LPMP {
   // standard visitor class for LPMP solver, when no custom visitor is given
   // do zrobienia: add xor arguments primalBoundComputationInterval, dualBoundComputationInterval with boundComputationInterval
   // do zrobienia: shall visitor depend on solver?
   //template<class SOLVER>
   class StandardVisitor {
      
      public:
      StandardVisitor(TCLAP::CmdLine& cmd)
         :
            posRealConstraint_(),
            posIntegerConstraint_(),
            maxIterArg_("","maxIter","maximum number of iterations of LPMP, default = 1000",false,1000,&posIntegerConstraint_,cmd),
            maxMemoryArg_("","maxMemory","maximum amount of memory (MB) LPMP is allowed to use",false,std::numeric_limits<INDEX>::max(),"positive integer",cmd),
            timeoutArg_("","timeout","time after which algorithm is stopped, in seconds, default = never, should this be type double?",false,std::numeric_limits<INDEX>::max(),&posIntegerConstraint_,cmd),
            // xor those //
            //boundComputationIntervalArg_("","boundComputationInterval","lower bound computation performed every x-th iteration, default = 5",false,5,"positive integer",cmd),
            primalComputationIntervalArg_("","primalComputationInterval","primal computation performed every x-th iteration, default = 5",false,5,&posIntegerConstraint_,cmd),
            primalComputationStartArg_("","primalComputationStart","iteration when to start primal computation, default = 1",false,1,&posIntegerConstraint_,cmd),
            lowerBoundComputationIntervalArg_("","lowerBoundComputationInterval","lower bound computation performed every x-th iteration, default = 1",false,1, &posIntegerConstraint_,cmd),
            ///////////////
            minDualImprovementArg_("","minDualImprovement","minimum dual improvement between iterations of LPMP",false,0.0,&posRealConstraint_,cmd),
            minDualImprovementIntervalArg_("","minDualImprovementInterval","the interval between which at least minimum dual improvement must occur",false,10,&posIntegerConstraint_,cmd),
            standardReparametrizationArg_("","standardReparametrization","mode of reparametrization",false,"anisotropic","{anisotropic|uniform}:${leave_percentage}",cmd),
            roundingReparametrizationArg_("","roundingReparametrization","mode of reparametrization for rounding primal solution:",false,"uniform:0.5","{anisotropic|uniform}:${leave_percentage}",cmd),
            primalTime_(0)
      {}

      template<typename LP_TYPE>
      LpControl begin(LP_TYPE& lp) // called, after problem is constructed. 
      {
         try {
            maxIter_ = maxIterArg_.getValue();
            maxMemory_ = maxMemoryArg_.getValue();
            remainingIter_ = maxIter_;
            minDualImprovement_ = minDualImprovementArg_.getValue();
            minDualImprovementInterval_ = minDualImprovementIntervalArg_.getValue();
            timeout_ = timeoutArg_.getValue();
            //boundComputationInterval_ = boundComputationIntervalArg_.getValue();
            primalComputationInterval_ = primalComputationIntervalArg_.getValue();
            primalComputationStart_ = primalComputationStartArg_.getValue();
            lowerBoundComputationInterval_ = lowerBoundComputationIntervalArg_.getValue();

            standard_reparametrization_ = lp_reparametrization(standardReparametrizationArg_.getValue());
            rounding_reparametrization_ = lp_reparametrization(roundingReparametrizationArg_.getValue());
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }

         beginTime_ = std::chrono::steady_clock::now();

         LpControl ret;
         ret.lp_repam = standard_reparametrization_;
         ret.computePrimal = false;
         ret.computeLowerBound = true;
         return ret;
      }

      // LpControl says what was last command to solver, return type gives next
      //template<typename SOLVER>
      LpControl visit(const LpControl c, const REAL lowerBound, const REAL primalBound)
      {
         lowerBound_.push_back(lowerBound); // rename to lowerBoundHistory_
         const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginTime_).count();

         // first output based on what lp solver did in last iteration
         if(c.computePrimal == false && c.computeLowerBound == false) {
            // output nothing
         } else { 
            if(verbosity >= 1) { 
              std::cout << "iteration = " << curIter_;
              if(c.computeLowerBound) {
                std::cout << ", lower bound = " << lowerBound;
              }
              if(c.computePrimal) {
                std::cout << ", upper bound = " << primalBound;
              }
              std::cout << ", time elapsed = " << timeElapsed/1000 << "." << (timeElapsed%1000)/10 << "s\n";
            }
         }

         curIter_++;
         remainingIter_--;

         LpControl ret;

         if(c.computePrimal) {
            prevLowerBound_ = curLowerBound_;
            curLowerBound_ = lowerBound;
         }
         // check if optimization has to be terminated
         if(remainingIter_ == 0) {
           if(verbosity >= 1) { std::cout << "One iteration remaining\n"; }
            ret.end = true;
            return ret;
         } 
         if(primalBound <= lowerBound + eps) {
            assert(primalBound + eps >= lowerBound);
            if(verbosity >= 1) { std::cout << "Primal cost " << primalBound << " greater equal lower bound " << lowerBound << "\n"; }
            ret.end = true;
            return ret;
         }
         if(timeout_ != std::numeric_limits<REAL>::max() && timeElapsed/1000 >= timeout_) {
            if(verbosity >= 1) { std::cout << "Timeout reached after " << timeElapsed << " seconds\n"; }
            remainingIter_ = std::min(INDEX(1),remainingIter_);
         }
         if(maxMemory_ > 0) {
            const INDEX memoryUsed = memory_used()/(1024*1024);
            if(maxMemory_ < memoryUsed) {
               remainingIter_ = std::min(INDEX(1),remainingIter_);
               if(verbosity >= 1) { std::cout << "Solver uses " << memoryUsed << " MB memory, aborting optimization\n"; }
            }
         }
         if(c.computeLowerBound && curIter_ >= minDualImprovementInterval_ && minDualImprovementArg_.isSet()) {
            assert(lowerBound_.size() >= minDualImprovementInterval_);
            const REAL prevLowerBound = lowerBound_[lowerBound_.size() - 1 - minDualImprovementInterval_];
            if(minDualImprovement_ > 0 && lowerBound - prevLowerBound < minDualImprovement_) {
               if(verbosity >= 1) { std::cout << "Dual improvement smaller than " << minDualImprovement_ << " after " << minDualImprovementInterval_ << " iterations, terminating optimization\n"; }
               remainingIter_ = std::min(INDEX(1),remainingIter_);
            }
         }

         if(remainingIter_ == 1) {
            ret.computePrimal = true;
            ret.computeLowerBound = true;
            ret.lp_repam = rounding_reparametrization_;
            return ret;
         }


         // determine next steps of solver
         ret.lp_repam = standard_reparametrization_;
         if(curIter_ >= primalComputationStart_ && (curIter_ - primalComputationStart_) % primalComputationInterval_ == 0) {
            ret.computePrimal = true; 
            ret.lp_repam = rounding_reparametrization_;
         }
         if(curIter_ % lowerBoundComputationInterval_ == 0) {
           ret.computeLowerBound = true;
         }
         return ret;
      }

      void end(const REAL lower_bound, const REAL upper_bound)
      {
         auto endTime = std::chrono::steady_clock::now();
         if(verbosity >= 1) { 
           std::cout << "final lower bound = " << lower_bound << ", upper bound = " << upper_bound << "\n";
           std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime_).count() << " milliseconds and " << curIter_ << " iterations.\n";
         }
      }
      
      using TimeType = decltype(std::chrono::steady_clock::now());
      TimeType GetBeginTime() const { return beginTime_; }
      //`REAL GetLowerBound() const { return curLowerBound_; }
      INDEX GetIter() const { return curIter_; }

      protected:
      PositiveRealConstraint posRealConstraint_;
      PositiveIntegerConstraint posIntegerConstraint_;
      // command line arguments TCLAP
      TCLAP::ValueArg<INDEX> maxIterArg_;
      TCLAP::ValueArg<INDEX> maxMemoryArg_;
      TCLAP::ValueArg<INDEX> timeoutArg_;
      //TCLAP::ValueArg<INDEX> boundComputationIntervalArg_;
      TCLAP::ValueArg<INDEX> primalComputationIntervalArg_;
      TCLAP::ValueArg<INDEX> primalComputationStartArg_;
      TCLAP::ValueArg<INDEX> lowerBoundComputationIntervalArg_;
      TCLAP::ValueArg<REAL> minDualImprovementArg_;
      TCLAP::ValueArg<INDEX> minDualImprovementIntervalArg_;
      TCLAP::ValueArg<std::string> standardReparametrizationArg_;
      TCLAP::ValueArg<std::string> roundingReparametrizationArg_;

      // command line arguments read out
      INDEX maxIter_;
      INDEX maxMemory_;
      INDEX timeout_;
      //INDEX boundComputationInterval_;
      INDEX primalComputationInterval_;
      INDEX primalComputationStart_;
      INDEX lowerBoundComputationInterval_;
      REAL minDualImprovement_;
      INDEX minDualImprovementInterval_;
      std::vector<REAL> lowerBound_; // do zrobienia: possibly make circular list out of this
      // do zrobienia: make enum for reparametrization mode
      lp_reparametrization standard_reparametrization_;
      lp_reparametrization rounding_reparametrization_;
      REAL rounding_reparametrization_leave_percentage_;

      // internal state of visitor
      INDEX remainingIter_;
      INDEX curIter_ = 0;
      REAL prevLowerBound_ = -std::numeric_limits<REAL>::max();
      REAL curLowerBound_ = -std::numeric_limits<REAL>::max();
      TimeType beginTime_;

      // primal
      //REAL bestPrimalCost_ = std::numeric_limits<REAL>::infinity();
      //REAL currentPrimalCost_ = std::numeric_limits<REAL>::infinity();

      //REAL bestLowerBound_ = -std::numeric_limits<REAL>::infinity();
      //REAL currentLowerBound_ = -std::numeric_limits<REAL>::infinity();

      //PrimalSolutionStorage currentPrimal_, bestPrimal_;
      INDEX primalTime_;
   };

   //template<class SOLVER>
   class StandardTighteningVisitor : public StandardVisitor //<SOLVER>
   {
      using BaseVisitorType = StandardVisitor; //<SOLVER>;
      public:
      StandardTighteningVisitor(TCLAP::CmdLine& cmd)
         :
            BaseVisitorType(cmd),
            tightenArg_("","tighten","enable tightening",cmd,false),
            tightenReparametrizationArg_("","tightenReparametrization","reparametrization mode used when tightening. Overrides primal computation reparametrization mode",false,"uniform:0.5","(uniform|anisotropic)",cmd),
            tightenIterationArg_("","tightenIteration","number of iterations after which tightening is performed for the first time, default = never",false,std::numeric_limits<INDEX>::max(),"positive integer", cmd),
            tightenIntervalArg_("","tightenInterval","number of iterations between tightenings",false,std::numeric_limits<INDEX>::max(),"positive integer", cmd),
            tightenConstraintsMaxArg_("","tightenConstraintsMax","maximal number of constraints to be added during tightening",false,20,"positive integer",cmd),
            tightenConstraintsPercentageArg_("","tightenConstraintsPercentage","maximal number of constraints to be added during tightening as percentage of number of initial factors",false,0.01,"positive real",cmd),
            posRealConstraint_(),
            tightenMinDualImprovementArg_("","tightenMinDualImprovement","minimum dual improvement after which to start tightening",false,std::numeric_limits<REAL>::infinity(),"positive real", cmd),
            tightenMinDualImprovementIntervalArg_("","tightenMinDualImprovementInterval","the interval between which at least minimum dual improvement may not occur for tightening",false,std::numeric_limits<INDEX>::max(), "positive integer", cmd),
            unitIntervalConstraint_(),
            tightenSlopeArg_("","tightenSlope","when slope of dual improvement becomes ${percentage} smaller than initial dual improvement slope, tighten", false, 1.0, &unitIntervalConstraint_, cmd)
            // do zrobienia: remove minDualIncrease and minDualDecreaseFactor
            //tightenMinDualIncreaseArg_("","tightenMinDualIncrease","obsolete: minimum increase which additional constraint must guarantee",false,0.0,&posRealConstraint_, cmd),
            //tightenMinDualDecreaseFactorArg_("","tightenMinDualDecreaseFactor","obsolete: factor by which to decrease minimum dual increase during tightening",false,0.5,&unitIntervalConstraint_, cmd)
      {
         //cmd.xorAdd(tightenConstraintsMaxArg_,tightenConstraintsPercentageArg_); // do zrobienia: this means that exactly one must be chosen. We want at most one to be chosen
      }

      template<typename LP_TYPE>
      LpControl begin(LP_TYPE& lp) // called, after problem is constructed. 
      {
         try {
            tighten_ = tightenArg_.getValue();
            tighten_reparametrization_ = lp_reparametrization(tightenReparametrizationArg_.getValue());
            tightenIteration_ = tightenIterationArg_.getValue();
            tightenInterval_ = tightenIntervalArg_.getValue();
            if(tightenConstraintsPercentageArg_.isSet() && tightenConstraintsMaxArg_.isSet()) {
               throw std::runtime_error("Only one of tightenConstraintsPercentage and tightenConstraintsMax may be set");
            }
            if(tightenConstraintsPercentageArg_.isSet()) {
               tightenConstraintsPercentage_ = tightenConstraintsPercentageArg_.getValue();
               tightenConstraintsMax_ = INDEX(tightenConstraintsPercentage_ * lp.number_of_factors());
            } else if(tightenConstraintsMaxArg_.isSet()) {
               tightenConstraintsMax_ = tightenConstraintsMaxArg_.getValue();
            } else if(tightenArg_.isSet()) {
               throw std::runtime_error("must set number of constraints to add");
            }
            tightenMinDualImprovement_ = tightenMinDualImprovementArg_.getValue();
            tightenMinDualImprovementInterval_ = tightenMinDualImprovementIntervalArg_.getValue();
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }

         return BaseVisitorType::begin(lp);
      }

      LpControl SetTighten(LpControl c)
      {
         c.tighten = true;
         c.tightenConstraints = tightenConstraintsMax_;
         c.lp_repam = tighten_reparametrization_;
         lastTightenIteration_ = this->GetIter();
         return c;
      }
      // the default
      //template<LPVisitorReturnType LP_STATE>
      LpControl visit(const LpControl c, const REAL lowerBound, const REAL primalBound)
      {
         auto ret = BaseVisitorType::visit(c, lowerBound, primalBound);

         if(tighten_) {
            iteration_after_tightening_++;
            const REAL cur_slope = std::max(lowerBound - prev_lower_bound_,REAL(0.0));
            if(iteration_after_tightening_ == 2) {
               tighten_slope_ = cur_slope;
            }
            if((this->GetIter() >= tightenIteration_ && 
                     (this->GetIter() >= lastTightenIteration_ + tightenInterval_ || 
                      (tightenSlopeArg_.isSet() && cur_slope < tightenSlopeArg_.getValue()*tighten_slope_)))) {
               if(verbosity >= 1) { std::cout << "Time to tighten\n"; }
               ret = SetTighten(ret);
               iteration_after_tightening_ = 0;
               tighten_slope_ = -std::numeric_limits<REAL>::infinity();
            } else if(this->GetIter() < tightenIteration_) {
               // check whether too small dual improvement necessitates tightening
               if(c.computeLowerBound && this->GetIter() > tightenMinDualImprovementInterval_ + lastTightenIteration_ && tightenMinDualImprovementArg_.isSet()) {
                  assert(this->lowerBound_.size() >= tightenMinDualImprovementInterval_);
                  const REAL prevLowerBound = lowerBound_[lowerBound_.size() - 1 - tightenMinDualImprovementInterval_];
                  if(tightenMinDualImprovement_ > 0 && lowerBound - prevLowerBound < tightenMinDualImprovement_) {
                     if(verbosity >= 1) { 
                       std::cout << "cur lower bound = " << lowerBound << " prev lower bound = " << prevLowerBound << "\n";
                       std::cout << "Dual improvement smaller than " << tightenMinDualImprovement_ << " after " << tightenMinDualImprovementInterval_ << " iterations, tighten\n";
                     }
                     ret = SetTighten(ret);
                     iteration_after_tightening_ = 0;
                     tighten_slope_ = -std::numeric_limits<REAL>::infinity();
                 }
              }
           }
         }
         prev_lower_bound_ = lowerBound;
         return ret;
      }
      /*
      INDEX Tighten()
      {
         auto tightenBeginTime = std::chrono::steady_clock::now();
         const INDEX constraintsAdded = BaseVisitorType::pd_.Tighten(tightenMinDualIncrease_, tightenConstraintsMax_);
         auto tightenEndTime = std::chrono::steady_clock::now();
         tightenTime_ += std::chrono::duration_cast<std::chrono::milliseconds>(tightenEndTime - tightenBeginTime).count();
         return constraintsAdded;
      }
      */

      protected:
      TCLAP::SwitchArg tightenArg_;
      TCLAP::ValueArg<std::string> tightenReparametrizationArg_;
      TCLAP::ValueArg<INDEX> tightenIterationArg_; // after how many iterations shall tightening be performed
      TCLAP::ValueArg<INDEX> tightenIntervalArg_; // interval between tightening operations.
      TCLAP::ValueArg<INDEX> tightenConstraintsMaxArg_; // How many constraints to add in tightening maximally
      TCLAP::ValueArg<REAL> tightenConstraintsPercentageArg_; // How many constraints to add in tightening maximally
      PositiveRealConstraint posRealConstraint_;
      PositiveIntegerConstraint posIntegerlConstraint_;
      TCLAP::ValueArg<REAL> tightenMinDualImprovementArg_;
      TCLAP::ValueArg<INDEX> tightenMinDualImprovementIntervalArg_;
      OpenUnitIntervalConstraint unitIntervalConstraint_;

      REAL prev_lower_bound_ = -std::numeric_limits<REAL>::infinity();
      TCLAP::ValueArg<REAL> tightenSlopeArg_;
      REAL tighten_slope_ = -std::numeric_limits<REAL>::infinity(); 
      INDEX iteration_after_tightening_ = 2; // this way tighten_slope will not be recomputed

      bool tighten_;
      lp_reparametrization tighten_reparametrization_;
      bool tightenInNextIteration_ = false;
      bool resumeInNextIteration_ = false;
      
      INDEX lastTightenIteration_ = 0;
      INDEX tightenIteration_;
      INDEX tightenInterval_;
      INDEX tightenConstraintsMax_;
      REAL tightenConstraintsPercentage_;
      REAL tightenMinDualImprovement_;
      INDEX tightenMinDualImprovementInterval_;
      REAL tightenMinDualIncrease_;
      REAL tightenMinDualDecreaseFactor_;
      
      INDEX tightenTime_ = 0;

   };
} // end namespace LPMP

#endif // LPMP_STANDARD_VISITOR_HXX
