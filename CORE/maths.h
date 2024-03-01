#ifndef MATHS_H
#define MATHS_H

#include <cmath>
#include <cstring>
#include <utility>
#include <vector>
#include <iostream>
#include <limits>

#include <array>

class Vec3d {
public:
  double v[3] = { };
  Vec3d() {}

  template <typename T>
  Vec3d(T &v1)
  {
    v[0]=v1[0]; v[1]=v1[1]; v[2]=v1[2];
  }
  Vec3d(const double x, const double y, const double z)
  {
    v[0]=x; v[1]=y; v[2]=z;
  }

  inline Vec3d operator-() const
  {
    return Vec3d(-v[0],-v[1],-v[2]);
  }
  inline double &operator[](size_t i)
  {
    return v[i];
  }
  inline double &operator()(size_t i)
  {
    return v[i];
  }
  inline const double &operator[](size_t i) const
  {
    return v[i];
  }
  inline const double &operator()(size_t i) const
  {
    return v[i];
  }
  template <typename T>
  inline Vec3d &operator+=(const T &v1)
  {
    v[0]+=v1[0];
    v[1]+=v1[1];
    v[2]+=v1[2];
    return *this;
  }
  template <typename T>
  inline Vec3d &operator-=(const T &v1)
  {
    v[0]-=v1[0];
    v[1]-=v1[1];
    v[2]-=v1[2];
    return *this;
  }
  inline Vec3d &operator*=(const double s)
  {
    v[0]*=s;
    v[1]*=s;
    v[2]*=s;
    return *this;
  }
  inline Vec3d operator-(const Vec3d &v2) const
  {
    return Vec3d((*this)[0]-v2[0], (*this)[1]-v2[1], (*this)[2]-v2[2]);
  }
  inline Vec3d operator+(const Vec3d &v2) const
  {
    return Vec3d((*this)[0]+v2[0], (*this)[1]+v2[1], (*this)[2]+v2[2]);
  }
  inline double operator*(const Vec3d &v2) const
  {
    return (*this)[0]*v2[0] + (*this)[1]*v2[1] + (*this)[2]*v2[2];
  }
  template <typename T>
  inline Vec3d operator*(const double s) const
  {
    return Vec3d((*this)[0]*s, (*this)[1]*s, (*this)[2]*s);
  }
  template <typename T>
  inline friend Vec3d operator*(const double s, const T &v)
  {
    return Vec3d(v[0]*s, v[1]*s, v[2]*s);
  }
  inline void set_zero() { v[0]=0; v[1]=0; v[2]=0; }

  inline friend std::ostream &operator<<(std::ostream &os, const Vec3d &v)
  {
    os << v[0] << "  " << v[1] << "  " << v[2];
    return os;
  }
  bool isApprox(const Vec3d &v1, double EPS) const{
    if (std::abs(v[0]-v1[0])>EPS) return false;
    if (std::abs(v[1]-v1[1])>EPS) return false;
    if (std::abs(v[2]-v1[2])>EPS) return false;
    return true;
  }
  inline bool operator==(const Vec3d &v1) const
  {
    double EPS = std::numeric_limits<double>::epsilon();
    return this->isApprox(v1, EPS);
  }
  inline bool operator!=(const Vec3d &v1) const
  {
    return !(*this==v1);
  }
  template <typename T>
  inline Vec3d &operator=(const T &v1)
  {
    v[0] = v1[0];
    v[1] = v1[1];
    v[2] = v1[2];
    return *this;
  }
  inline void operator()(const double x, const double y, const double z)
  {
    v[0]=x;
    v[1]=y;
    v[2]=z;
  }
  inline Vec3d cross(const Vec3d &v) const
  {
    Vec3d temp;
    temp[0] = (*this)[1] * v[2] - (*this)[2] * v[1];
    temp[1] = -((*this)[0] * v[2] - (*this)[2] * v[0]);
    temp[2] = (*this)[0] * v[1] - (*this)[1] * v[0];
    return temp;
  }

  inline double norm() const {
    return std::sqrt((*this)*(*this));
  }
};

class Mat6R3C {
  // !! Row major order !!
public:
  int rows=6;
  int cols=3;
  double m[18];
  inline double &operator()(const size_t i, const size_t j)
  {
    return m[i*cols + j];
  }
  inline Vec3d row(const size_t r) {
    return Vec3d(m[r*cols],m[r*cols+1],m[r*cols+2]);
  }
};



class aed_type2 {
  // Simple vector class
  //
  // Operator * is only defined for the dot product
  // so either v1*v2 and v2*v1 always give dot product
  size_t size_= 0;
  size_t stride_= 1;
  bool allocated_=false;  // false is required for viewing
public:
  double *data_ = NULL;
  aed_type2() {}
  aed_type2(size_t n, size_t s=1):
    size_(n),
    stride_(s),
    allocated_(true),
    data_((double*)malloc(s*n*sizeof(data_)))
  {}
  ~aed_type2() {
    if (allocated_)
      free(data_);
  }
  aed_type2(const aed_type2 &other) {
    // the other might be a view vector

    size_= other.size_;
    stride_ = other.stride_;

    //if (size_ > 0) {
      data_ = (double*)malloc(stride_*size_*sizeof(data_));
      std::memcpy(data_, other.data_, stride_*size_*sizeof(data_));
      allocated_ = true;
    //}
  }
  aed_type2 &operator=(aed_type2 other) {
    std::swap(size_, other.size_);
    std::swap(allocated_, other.allocated_);
    std::swap(stride_, other.stride_);
    std::swap(data_, other.data_);
    return *this;
  }
  template <typename T>
  aed_type2(T &v)
  {
    resize(v.size());
    size_ = v.size();
    //TODO stride???
    for (size_t i=0; i<v.size(); ++i) {
      data_[i]=v[i];
    }
  }

  /** This constructor does not allocate memory */
  aed_type2(double* arr, size_t size, size_t stride=1);

  /** Allocated if copy is true otherwise provides a view to arr */
  aed_type2(double* arr, size_t size, bool copy, size_t stride=1);

  //aed_type2(double* arr, size_t size, bool copy=true, size_t stride=1);
  inline double &operator[](size_t i)
  {
    return data_[i*stride_];
  }
  inline double &operator()(size_t i)
  {
    return data_[i*stride_];
  }
  inline const double &operator[](size_t i) const
  {
    return data_[i*stride_];
  }
  inline const double &operator()(size_t i) const
  {
    return data_[i*stride_];
  }
  //inline aed_type2 transpose() const
  //{
  //    return *this;
  //}
  // Maths
  //inline aed_type2 operator-() const
  //{
  //  aed_type2 temp(*this);
  //  for (size_t i=0;i<size();++i) {
  //    temp[i]*=-1;
  //  }
  //  return temp;;
  //}
  template <typename T>
  inline aed_type2 &operator+=(const T &v1)
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]+=v1[i];
    }
    return *this;
  }
  inline aed_type2 &operator+=(const double d)
  // element wise increment
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]+=d;
    }
    return *this;
  }
  inline aed_type2 operator+(const double d) const
  // element wise increment
  {
    aed_type2 temp(*this);
    for (size_t i=0;i<size();++i) {
      temp[i]+=d;
    }
    return temp;
  }
  template <typename T>
  inline aed_type2 &operator-=(const T &v1)
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]-=v1[i];
    }
    return *this;
  }
  inline aed_type2 &operator-=(const double d)
  // element wise decrement
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]-=d;
    }
    return *this;
  }
  inline aed_type2 operator-(const double d) const
  // element wise decrement
  {
    aed_type2 temp(*this);
    for (size_t i=0;i<size();++i) {
      temp[i]-=d;
    }
    return temp;
  }
  inline aed_type2 &operator*=(const double s)
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]*=s;
    }
    return *this;
  }
  template <typename T>
  inline aed_type2 operator-(const T &v2) const
  {
    aed_type2 temp(*(this));
    for (size_t i=0;i<size();++i) {
      temp[i]-=v2[i];
    }
    return temp;
  }
  template <typename T>
  inline aed_type2 operator+(const T &v2) const
  {
    aed_type2 temp((*this));
    for (size_t i=0;i<size();++i) {
      temp[i]+=v2[i];
    }
    return temp;
  }
  template <typename T>
  inline double operator*(const T &v2) const
  {
    double temp=0;
    for (size_t i=0;i<size();++i) {
      temp+=(*this)(i)*v2(i);
    }
    return temp;
  }
  inline aed_type2 operator*(const double s) const
  {
    aed_type2 temp(*(this));
    for (size_t i=0;i<size();++i) {
      temp[i]*=s;
    }
    return temp;
  }
  inline friend aed_type2 operator*(const double s,
                                    const aed_type2 &v)
  {
    aed_type2 temp(v);
    for (size_t i=0;i<v.size();++i) {
      temp[i]*=s;
    }
    return temp;
  }
  inline void set_zero()
  {
    for (size_t i=0;i<size();++i) {
      (*this)[i]=0;
    }
  }

  size_t size() const
  {
    return size_;
  }
  size_t rows() const
  {
    return size_;
  }
  bool isApprox(const aed_type2 &v1, double EPS) const{
    for (size_t i=0;i<size();++i) {
      if (std::abs((*this)[i]-v1[i])>EPS) return false;
    }
    return true;
  }
  inline bool operator==(const aed_type2 &v1) const
  {
    double EPS = std::numeric_limits<double>::epsilon();
    return this->isApprox(v1, EPS);
  }
  inline bool operator!=(const aed_type2 &v1) const
  {
    return !(*this==v1);
  }
  inline void resize(const size_t n) {
    if (n != size_) {
      if (allocated_) {
        data_ = (double*)realloc(data_, stride_*n*sizeof(data_));
        size_ = n;
      }
      else {
        data_ = (double*)malloc(stride_*n*sizeof(data_));
        size_ = n;
        allocated_ = true;
      }
    }
  }
  double *data() {
    return data_;
  }
  const double *data() const {
    return data_;
  }
  aed_type2 copy(size_t nelements, size_t start=0) {
    if (nelements+start>size())
      throw std::runtime_error("Out of range\n");
    aed_type2 temp(nelements);
    for (size_t i=0; i<nelements; ++i) {
      temp[i]=(*this)[start+i];
    }
    return temp;
  }
  /** Return the pointer to the first element of the data array */
  inline double *ptr() {
    return &data()[0];
  }
  inline double mean() const {
    double mean=0;
    for (size_t i=0; i<size(); ++i)
      mean+=(*this)(i);
    return mean/size();
  }
  inline double std_dev(double mean, double N) const {
    double stdev=0;
    for (size_t i=0; i<size(); ++i)
      stdev += (((*this)(i)-mean))*((*this)(i)-mean);
    return std::sqrt(stdev/N);
  }
};

inline std::ostream &operator<<(std::ostream &os, const aed_type2 &v)
{
  for (size_t i=0;i<v.size()-1;++i) {
    os << v[i] << "  ";
  }
  os << v[v.size()-1];
  return os;
}

class aeds_type2 {
  // column-wise
public:
  std::vector<aed_type2> data_;
  aeds_type2() {}
  aeds_type2(size_t dim, size_t natoms):
    data_(natoms, aed_type2(dim))
  {}
  inline double &operator()(size_t i, size_t j)
  {
    return data_[i][j];
  }
  inline aed_type2 &operator[](size_t i)
  {
    return data_[i];
  }
  size_t cols() const
  {
    return data_.size();
  }
  size_t rows() const
  {
    return data_[0].size();
  }
  size_t size() const
  {
    return data_.size();
  }
  inline aed_type2 &col(size_t i)
  {
    return data_[i];
  }
  inline const aed_type2 &col(size_t i) const
  {
    return data_[i];
  }
  inline void resize(size_t dim, size_t natoms)
  {
    data_.resize(natoms);
    for (size_t i = 0; i < natoms; ++i)
      data_[i].resize(dim);
  }
  inline void set_zero()
  {
    for (auto &a:data_) a.set_zero();
  }
  inline friend std::ostream &operator<<(std::ostream &os, const aeds_type2 &v)
  {
    for (const aed_type2 &a:v.data_) os << a << std::endl;
    return os;
  }
};

class fd_type {
  // column-wise

  size_t rows_=0;
  std::vector<double> data;
  aed_type2 aed_type2_ptr[3];
public:
  fd_type();
  fd_type(size_t dim);
  inline double &operator()(size_t i, size_t k)
  {
    return data[3*i+k];
  }
  inline const double &operator()(size_t i, size_t k) const
  {
    return data[3*i+k];
  }
  inline aed_type2 &operator[](size_t k)
  {
    return aed_type2_ptr[k];
  }
  inline const aed_type2 &operator()(size_t k) const
  {
    return aed_type2_ptr[k];
  }
  size_t rows() const
  {
    return rows_;
  }
  inline void set_zero()
  {
    for (size_t i=0; i<3*rows_; ++i) data[i]=0;
  }
  inline friend std::ostream &operator<<(std::ostream &os, const fd_type &v)
  {
    os << v.aed_type2_ptr[0] << std::endl;
    os << v.aed_type2_ptr[1] << std::endl;
    os << v.aed_type2_ptr[2] << std::endl;
    return os;
  }
  template <typename T>
  inline Vec3d operator*(const T &v2) const
  {
    Vec3d temp;
    for (size_t k=0; k<3; ++k) {
      for (size_t i=0;i<v2.size();++i) {
        temp[k]+=(*this)(i,k)*v2[i];
      }
    }
    return temp;
  }
  inline fd_type operator-(const fd_type &v2) const
  {
    fd_type temp(v2.rows());
    for (size_t i=0;i<3*rows();++i)
      temp.data[i] = data[i]-v2.data[i];
    return temp;
  }
  inline void operator*=(double s) 
  {
    for (size_t i=0;i<3*rows();++i)
      data[i]*=s;
  }
  //inline fd_type operator-() const
  //{
  //  fd_type temp(*this);
  //  for (size_t k=0; k<3; ++k) {
  //    for (size_t i=0;i<temp(k).size();++i) {
  //      temp(i,k)*=-1;
  //    }
  //  }
  //  return temp;
  //}
};



template <typename D>
class Matrix_Base {
  // column-major order
  size_t rows_;
  size_t cols_;

public:
  D m_;
  Matrix_Base<D>() {};
  Matrix_Base<D>(D &d, size_t r, size_t c):
    rows_(r),
    cols_(c),
    m_(d)
  {}
  Matrix_Base<D>(double *d, size_t r, size_t c):
    rows_(r),
    cols_(c),
    m_(d,d+r*c)
  {}
  Matrix_Base<D>(size_t r, size_t c):
    rows_(r),
    cols_(c),
    m_(r*c)
  {}

  inline double& operator()(size_t i, size_t j)
  {
    return m_[i + rows_*j];
  }

  inline double operator()(size_t i, size_t j) const
  {
    return m_[i + rows_*j];
  }
  inline double& operator[](size_t i)
  {
    return m_[i];
  }

  inline double operator[](size_t i) const
  {
    return m_[i];
  }

  inline size_t rows() const {
    return rows_;
  }
  inline size_t cols() const {
    return cols_;
  }
  inline double *data() {
    return m_.data();
  }
  const double *data() const {
    return m_.data();
  }
  /** Return the pointer to the first element of the data array */
  inline double *ptr() {
    return &data()[0];
  }
  inline bool exists(size_t row, size_t col) const
  {
    return (row < rows_ && col < cols_);
  }
  /** Provide a view of the r row of the matrix */
  inline aed_type2 row(size_t r) {
    size_t stride = rows();
    aed_type2 temp(&m_.data()[r], cols(), stride);
    return temp;
  }
  /** Provide a view of the c column of the matrix */
  inline aed_type2 col(size_t c) {
    size_t start = rows()*c;
    aed_type2 temp(&m_.data()[start], rows());
    return temp;
  }
  inline void resize(size_t r, size_t c)
  {
    m_.resize(r*c);
    rows_= r;
    cols_= c;
  }
  inline void set_zero() {
    for (size_t i=0;i<m_.size();++i) {
      m_[i]=0;
    }
  }
  inline size_t size() const {
    // Number of coefficients in the matrix
    return m_.size();
  }
  //inline t_pred operator*(const aed_type2) {
  //}
  inline aed_type2 diagonal() const {
    aed_type2 temp(rows());
    for (size_t i=0; i<rows(); ++i) {
      temp(i)=(*this)(i,i);
    }
    return temp;
  }
  inline friend std::ostream &operator<<(std::ostream &os,
                                         const Matrix_Base<D> &m)
  {
    for (size_t r=0;r<m.rows();++r) {
      for (size_t c=0;c<m.cols()-1;++c) {
        os << m(r,c) << "  ";
      }
      os << m(r,m.cols()-1) << std::endl;
    }
    return os;
  }
  template <typename T>
  inline aed_type2 operator*(const T &v) const
  {
    // naive matrix vector product
    //m is M x N matrix and v is N x 1 column vector
    aed_type2 temp((*this).rows());
    temp.set_zero();
    for (size_t i=0; i<(*this).rows(); i++){
      for (size_t j=0; j<(*this).cols(); j++) {
        temp(i) += (*this)(i,j)*v(j);
      }
    }
    return temp;
  }
  inline void set_rows(const size_t r)
  {
    rows_=r;
  }
  inline void set_cols(const size_t c)
  {
    cols_=c;
  }
  // Maths
  inline Matrix_Base<D> operator-() const
  {
    Matrix_Base<D> temp(*this);
    for (size_t i=0;i<temp.size();++i) {
      temp.data()[i]*=-1;
    }
    return temp;;
  }
  template <typename T>
  inline Matrix_Base<D> &operator+=(const T &m)
  {
    for (size_t i=0;i<m_.size();++i) {
      m_[i]+=m.data()[i];
    }
    return *this;
  }
  template <typename T>
  inline Matrix_Base<D> &operator-=(const T &m)
  {
    for (size_t i=0;i<m_.size();++i) {
      m_[i]-=m.data()[i];
    }
    return *this;
  }
  inline double trace() const {
    double temp=0;
    for (size_t i=0; i<rows(); ++i)
      temp+= (*this)(i,i);
    return temp;
  }
  bool isApprox(const Matrix_Base<D> &m, double EPS) const {
    for (size_t i=0;i<size();++i) {
      if (std::abs(data()[i]-m.data()[i])>EPS) return false;
    }
    return true;
  }
  inline bool operator==(const aed_type2 &v1) const
  {
    double EPS = std::numeric_limits<double>::epsilon();
    return this->isApprox(v1, EPS);
  }
  inline aed_type2 col(size_t c) const {
    size_t start = rows()*c;

    aed_type2 temp(rows());
    for (size_t i=0; i<rows(); ++i) {
      temp(i)=data()[i+start];
    }
    return temp;
  }
};
template<typename T=std::vector<double>>
class Matrix_T: public Matrix_Base<T> {
public:
  Matrix_T() {};
  Matrix_T(T &d, size_t r, size_t c):
    Matrix_Base<T>(d,r,c)
  {}
  Matrix_T(double *d, size_t r, size_t c):
    Matrix_Base<T>(d,r,c)
  {};
  Matrix_T(size_t r, size_t c):
    Matrix_Base<T>(r,c)
  {}
};
typedef Matrix_T<> Matrix;

template<typename T=std::array<double,9>>
class Matrix3d_T: public Matrix_Base<T> {
  // class for 3x3 array
public:
  Matrix3d_T()
  {
    Matrix_Base<T>::set_rows(3);
    Matrix_Base<T>::set_cols(3);
  }
  Matrix3d_T(T &d, size_t r, size_t c):
    Matrix_Base<T>(d,r,c)
  {}
  inline Vec3d row(size_t r) const
  {
    Vec3d temp((*this)[r], (*this)[r+3], (*this)[r+6]);
    return temp;
  }
  inline Vec3d col(size_t c) const
  {
    size_t s=3*c;
    Vec3d temp((*this)[s], (*this)[s+1], (*this)[s+2]);
    return temp;
  }
  inline Matrix3d_T<T> transpose() {
    Matrix3d_T<T> temp;
    for (size_t i=0; i<3; i++){
      for (size_t j=0; j<3; j++) {
        temp(j,i) = (*this)(i,j);
      }
    }
    return temp;
  }
  template <typename D>
  inline Vec3d operator*(const D &v) const
  {
    // naive matrix vector product
    //m is M x N matrix and v is N x 1 column vector
    Vec3d temp;
    temp.set_zero();
    for (size_t i=0; i<3; i++){
      for (size_t j=0; j<3; j++) {
        temp(i) += (*this)(i,j)*v(j);
      }
    }
    return temp;
  }
  inline Matrix3d_T<T> inverse() const
  {
    // computes the inverse of a 3x3 matrix
    double det = (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) -
      (*this)(2, 1) * (*this)(1, 2)) -
      (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) -
      (*this)(1, 2) * (*this)(2, 0)) +
      (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) -
      (*this)(1, 1) * (*this)(2, 0));

    double invdet = 1 / det;

    Matrix3d_T<T> minv;
    minv(0, 0) = ((*this)(1, 1) * (*this)(2, 2) -
      (*this)(2, 1) * (*this)(1, 2)) * invdet;
    minv(0, 1) = ((*this)(0, 2) * (*this)(2, 1) -
      (*this)(0, 1) * (*this)(2, 2)) * invdet;
    minv(0, 2) = ((*this)(0, 1) * (*this)(1, 2) -
      (*this)(0, 2) * (*this)(1, 1)) * invdet;
    minv(1, 0) = ((*this)(1, 2) * (*this)(2, 0) -
      (*this)(1, 0) * (*this)(2, 2)) * invdet;
    minv(1, 1) = ((*this)(0, 0) * (*this)(2, 2) -
      (*this)(0, 2) * (*this)(2, 0)) * invdet;
    minv(1, 2) = ((*this)(1, 0) * (*this)(0, 2) -
      (*this)(0, 0) * (*this)(1, 2)) * invdet;
    minv(2, 0) = ((*this)(1, 0) * (*this)(2, 1) -
      (*this)(2, 0) * (*this)(1, 1)) * invdet;
    minv(2, 1) = ((*this)(2, 0) * (*this)(0, 1) -
      (*this)(0, 0) * (*this)(2, 1)) * invdet;
    minv(2, 2) = ((*this)(0, 0) * (*this)(1, 1) -
      (*this)(1, 0) * (*this)(0, 1)) * invdet;

    return minv;
  }
  template <typename D>
  inline Matrix3d_T<T> &operator+=(const D &m)
  {
    for (size_t i=0;i<9;++i) {
      Matrix_Base<T>::m_[i]+=m.data()[i];
    }
    return *this;
  }
  template <typename D>
  inline Matrix3d_T<T> &operator-=(const D &m)
  {
    for (size_t i=0;i<9;++i) {
      Matrix_Base<T>::m_[i]-=m.data()[i];
    }
    return *this;
  }
};

typedef Matrix3d_T<> Matrix3d;



#endif
