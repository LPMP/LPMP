#ifndef LPMP_THREE_DIMENSIONAL_VARIABLE_ARRAY_HXX
#define LPMP_THREE_DIMENSIONAL_VARIABLE_ARRAY_HXX

#include <vector>
#include <array>
#include <cassert>

namespace LPMP {

// given fixed first index, the matrices can be of differing sizes
template<typename T>
class three_dimensional_variable_array {
public:

    three_dimensional_variable_array() {}

    template<typename ITERATOR>
    std::size_t set_dimensions(ITERATOR size_begin, ITERATOR size_end)
    {
        dimensions_.clear();
        dimensions_.reserve(std::distance(size_begin, size_end));
        std::size_t offset = 0;
        std::size_t size = 0;
        for(std::size_t i=0; i<std::distance(size_begin, size_end); ++i) {
           dimensions_.push_back( dim_entry(size_begin[i][0], size_begin[i][1], offset) );
           offset += dimensions_[i][0] * dimensions_[i][1];
           size += dimensions_[i][0] * dimensions_[i][1];
        } 
        return size;
    }

    template<typename ITERATOR>
    three_dimensional_variable_array(ITERATOR size_begin, ITERATOR size_end)
    {
       const auto s = set_dimensions(size_begin, size_end);
       data_.resize(s);
    }

    template<typename ITERATOR>
    three_dimensional_variable_array(ITERATOR size_begin, ITERATOR size_end, T value)
    {
       const auto s = set_dimensions(size_begin, size_end);
       data_.resize(s, value);
    }

    template<typename ITERATOR>
    void resize(ITERATOR size_begin, ITERATOR size_end, T value)
    {
       const auto s = set_dimensions(size_begin, size_end);
       data_.resize(s, value);
    }

    template<typename ITERATOR>
    void resize(ITERATOR size_begin, ITERATOR size_end)
    {
       const auto s = set_dimensions(size_begin, size_end);
       data_.resize(s);
    }

    T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) { 
        assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
        const auto offset = dimensions_[x1].offset;
        return data_[offset + x2*dim3(x1) + x3]; 
    }
    const T& operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const 
    {
        assert(x1<dim1() && x2<dim2(x1) && x3<dim3(x1));
        const auto offset = dimensions_[x1].offset;
        return data_[offset + x2*dim3(x1) + x3]; 
    }      

    const std::size_t dim1() const { return dimensions_.size(); }
    const std::size_t size() const { return dim1(); }
    const std::size_t dim2(const std::size_t x1) const
    {
        assert(x1<dim1());
        return dimensions_[x1][0];
    }
    const std::size_t dim3(const std::size_t x1) const
    {
        assert(x1<dim1());
        return dimensions_[x1][1];
    }

    std::vector<T>& data() { return data_; }
    const std::vector<T>& data() const { return data_; }

    struct const_matrix_access_object {
        const_matrix_access_object(const T* data, const std::size_t dim1, const std::size_t dim2)
        : data_(data),
        dim1_(dim1),
        dim2_(dim2)
        {}

        const std::size_t dim1() const { return dim1_; }
        const std::size_t dim2() const { return dim2_; }
        const std::size_t size() const { return dim1()*dim2(); }
        const T& operator()(const std::size_t x1, const std::size_t x2) const 
        { 
            assert(x1<dim1() && x2<dim2());
            return data_[x1*dim2() + x2]; 
        }

        private:
        const T* data_;
        const std::size_t dim1_;
        const std::size_t dim2_; 
    };
    struct matrix_access_object {
        matrix_access_object(T* data, const std::size_t dim1, const std::size_t dim2)
        : data_(data),
        dim1_(dim1),
        dim2_(dim2)
        {}

        const std::size_t dim1() const { return dim1_; }
        const std::size_t dim2() const { return dim2_; }
        const std::size_t size() const { return dim1()*dim2(); }
        const T& operator()(const std::size_t x1, const std::size_t x2) const 
        { 
            assert(x1<dim1() && x2<dim2());
            return data_[x1*dim2() + x2]; 
        }
        T& operator()(const std::size_t x1, const std::size_t x2) 
        {
            assert(x1<dim1() && x2<dim2());
            return data_[x1*dim2() + x2]; 
        }

        private:
        T* data_;
        const std::size_t dim1_;
        const std::size_t dim2_; 
    };

    const_matrix_access_object operator[](const std::size_t x1) const
    {
        assert(x1<dim1());
        const auto offset = dimensions_[x1].offset;
        return const_matrix_access_object(&(data_[offset]), dimensions_[x1][0], dimensions_[x1][1]);
    }
    matrix_access_object operator[](const std::size_t x1)
    {
        assert(x1<dim1());
        const auto offset = dimensions_[x1].offset;
        return matrix_access_object(&(data_[offset]), dimensions_[x1][0], dimensions_[x1][1]);
    }

    auto size_begin() const { return dimensions_.begin(); }
    auto size_end() const { return dimensions_.end(); }

private:
    template<typename ITERATOR>
    static std::size_t size(ITERATOR size_begin, ITERATOR size_end)
    {
        std::size_t s=0;
        for(auto size_it=size_begin; size_it!=size_end; ++size_it) {
            s += (*size_it)[0] * (*size_it)[1];
        }
        return s;
    }
    struct dim_entry: public std::array<std::size_t,2> { // dim2, dim3
       dim_entry(const std::size_t dim_1, const std::size_t dim_2, const std::size_t _offset) : std::array<std::size_t,2>({dim_1,dim_2}), offset(_offset) {}
       std::size_t offset; // offset into data array
    };

    std::vector<dim_entry> dimensions_; 
    std::vector<T> data_; 
};

} // namespace LPMP

#endif // LPMP_THREE_DIMENSIONAL_VARIABLE_ARRAY_HXX

