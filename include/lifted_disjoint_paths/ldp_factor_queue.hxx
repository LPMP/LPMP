#ifndef LDP_FACTOR_QUEUE_HXX
#define LDP_FACTOR_QUEUE_HXX

#include <queue>
#include <cmath>
#include <config.hxx>


namespace LPMP {



template <class T>
class ForQueueComparator{
public:
    bool operator()(const std::tuple<double,size_t,T*>& lhs, const  std::tuple<double,size_t,T*>& rhs) {
        if(std::abs(std::get<0>(rhs)-std::get<0>(lhs))<eps){
            return std::get<1>(rhs) < std::get<1>(lhs);
        }
        else{
            return std::get<0>(lhs)<std::get<0>(rhs);
        }
    }
};


template<class T>
class LdpFactorQueue{
public:

    using QueueElement=std::tuple<double,size_t,T*>;

    LdpFactorQueue(){
        numberOfElements=0;
    }

    void clearPriorityQueue();
    void insertToQueue(const double& improvementValue,T* pFactor );

    const T* getTopFactorPointer() const {
        assert(!stableQueue.empty());
        return std::get<2>(stableQueue.top());

    }

    const double& getTopImprovementValue() const {
        assert(!stableQueue.empty());
        return std::get<0>(stableQueue.top());
    }

    bool isQueueEmpty() const {
        return stableQueue.empty();
    }

    void removeTopElement(){
        assert(!stableQueue.empty());
        T* p=std::get<2>(stableQueue.top());
        delete p;
        stableQueue.pop();
        numberOfElements--;
    }

private:

std::priority_queue<QueueElement,std::vector<QueueElement>,ForQueueComparator<T>> stableQueue;
size_t numberOfElements;

};



template <class T>
inline void LdpFactorQueue<T>::insertToQueue(const double& improvementValue,T* pFactor){
    stableQueue.push(std::make_tuple(improvementValue,numberOfElements,pFactor));
    numberOfElements++;
}



template <class T>
inline void LdpFactorQueue<T>::clearPriorityQueue() {
    while(!stableQueue.empty()){
       removeTopElement();
    }
    assert(numberOfElements==0);
}



}
#endif // LDP_FACTOR_QUEUE_HXX
