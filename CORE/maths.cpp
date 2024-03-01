#include "maths.h"
#include "lapack.h"


// aed_type2 class

aed_type2::aed_type2(double* arr, size_t size, size_t stride)
 {
     stride_ = stride;
     size_   = size;
     data_   = arr;
 }

aed_type2::aed_type2(double* arr, size_t size, bool copy, size_t stride)
{
    stride_= stride;
    size_ = size;
    if (copy) {
        data_ = (double*)malloc(stride_*size_*sizeof(data_));
        std::memcpy(data_, arr, stride_*size_*sizeof(data_));
    }
    else {
        data_ = arr;
    }
}

// fd_type
fd_type::fd_type() 
{}

fd_type::fd_type(size_t dim):
  rows_(dim),
  data(3*dim)
{
  aed_type2_ptr[0] = aed_type2(&data[0],dim,(size_t)1);
  aed_type2_ptr[1] = aed_type2(&data[dim],dim,(size_t)1);
  aed_type2_ptr[2] = aed_type2(&data[2*dim],dim,(size_t)1);
}

