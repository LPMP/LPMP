#ifndef LPMP_SPINLOCK_HXX
#define LPMP_SPINLOCK_HXX 

#include <atomic>

#if defined(_MSC_VER) && _MSC_VER >= 1310 && ( defined(_M_IX86) || defined(_M_X64) )

extern "C" void _mm_pause();

#define LPMP_SMT_PAUSE _mm_pause();

#elif defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )

#define LPMP_PAUSE __asm__ __volatile__( "rep; nop" : : : "memory" );

#endif


namespace LPMP { 

  class spinlock {
    std::atomic_flag locked = ATOMIC_FLAG_INIT ;
    public:
    void lock() {
      while (locked.test_and_set(std::memory_order_acquire)) { 
        LPMP_PAUSE
      }
  }
  void unlock() {
    locked.clear(std::memory_order_release);
  }
};

} // end namespace LPMP

#endif // LPMP_SPINLOCK_HXX
