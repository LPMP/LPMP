#pragma once

template<class T>
    class CopyableAtomic : public std::atomic<T>
    {
    public:
        CopyableAtomic() :
            std::atomic<T>(T{})
        {}

        constexpr CopyableAtomic(T desired): std::atomic<T>(desired){}

        constexpr CopyableAtomic(const CopyableAtomic<T>& other): CopyableAtomic(other.load(std::memory_order_relaxed))
        {}

        CopyableAtomic& operator=(const CopyableAtomic<T>& other) {
            this->store(other, std::memory_order_relaxed);
            return *this;
        }

        CopyableAtomic& operator=(const T& other) {
            this->store(other, std::memory_order_relaxed);
            return *this;
        }
    };