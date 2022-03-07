#ifndef LPMP_SPINLOCK_HXX
#define LPMP_SPINLOCK_HXX 

#include <atomic>

namespace LPMP {

class spinlock {
    private:
        std::atomic<bool> lock_ = {false};

    public:
        void lock() { while(lock_.exchange(true, std::memory_order_acquire)); }

        void unlock() { lock_.store(false, std::memory_order_release); }
};

} // end namespace LPMP

#endif // LPMP_SPINLOCK_HXX
