#include "complex_arith.h"

// add complex numbers
HOSTDEVICE cudacomplex operator+(const cudacomplex REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a.value.x + b.value.x, a.value.y  + b.value.y }};
   return result;
}

/*
// increase complex number by another complex number
HOSTDEVICE cudacomplex operator+=(const cudacomplex REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a.value.x + b.value.x, a.value.y  + b.value.y }};
   return result;
}

// increase complex number by a scalar
HOSTDEVICE cudacomplex operator+=(const cudacomplex REF(a), const flt_doub REF(b)) {
   cudacomplex result = {{ a.value.x + b, a.value.y }};
   return result;
}
*/
// add scalar to complex
HOSTDEVICE cudacomplex operator+(const cudacomplex REF(a), const flt_doub REF(b)) {
   cudacomplex result = {{ a.value.x + b, a.value.y }};
   return result;
}

// add complex to scalar
HOSTDEVICE cudacomplex operator+(const flt_doub REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a + b.value.x, b.value.y }};
   return result;
}

// subtract complex numbers
HOSTDEVICE cudacomplex operator-(const cudacomplex REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a.value.x - b.value.x, a.value.y  - b.value.y }};
   return result;
}

// negate a complex number
HOSTDEVICE cudacomplex operator-(const cudacomplex REF(a)) {
   cudacomplex result = {{ -a.value.x, -a.value.y }};
   return result;
}

// subtract scalar from complex
HOSTDEVICE cudacomplex operator-(const cudacomplex REF(a), const flt_doub REF(b)) {
   cudacomplex result = {{ a.value.x - b, a.value.y }};
   return result;
}

// subtract complex from scalar
HOSTDEVICE cudacomplex operator-(const flt_doub REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a - b.value.x, -b.value.y }};
   return result;
}

// multiply complex numbers
HOSTDEVICE cudacomplex operator*(const cudacomplex REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a.value.x * b.value.x - a.value.y * b.value.y,
                           a.value.y * b.value.x + a.value.x * b.value.y }};
   return result;
}

/*
DEVICE cudacomplex cexp_dev(const cudacomplex REF(a)) {
      cudacomplex result = {{__powf(M_E,a.value.x)*cosf(a.value.y), __powf(M_E,a.value.x)*sinf(a.value.y)}};
      return result;
}
*/

HOSTDEVICE cudacomplex cexp(const cudacomplex REF(a)) {
      cudacomplex result = {{powf(M_E,a.value.x)*cosf(a.value.y), powf(M_E,a.value.x)*sinf(a.value.y)}};
      return result;
}

HOSTDEVICE double carg(const cudacomplex REF(a)) {
      double result = atan2(a.value.y, a.value.x); //{{powf(M_E,a.value.x)*cosf(a.value.y), powf(M_E,a.value.x)*sinf(a.value.y)}};
      return result;
}

DEVICE cudacomplex cexp_dev(const cudacomplex REF(a)) {
      cudacomplex result = {{__powf(M_E,a.value.x)*__cosf(a.value.y), __powf(M_E,a.value.x)*__sinf(a.value.y)}};
      return result;
}
// multiply complex with scalar
HOSTDEVICE cudacomplex operator*(const cudacomplex REF(a), const flt_doub REF(b)) {
   cudacomplex result = {{ a.value.x * b, a.value.y * b }};
   return result;
}

// multiply scalar with complex
HOSTDEVICE cudacomplex operator*(const flt_doub REF(a), const cudacomplex REF(b)) {
   cudacomplex result = {{ a * b.value.x, a * b.value.y }};
   return result;
}

// divide complex numbers
HOSTDEVICE cudacomplex operator/(const cudacomplex REF(a), const cudacomplex REF(b)) {
   flt_doub tmp = ( b.value.x * b.value.x + b.value.y * b.value.y );
   cudacomplex result = {{ (a.value.x * b.value.x + a.value.y * b.value.y ) / tmp,
                           (a.value.y * b.value.x - a.value.x * b.value.y ) / tmp }};
   return result;
}
// divide complex by scalar
HOSTDEVICE cudacomplex operator/(const cudacomplex REF(a), const flt_doub REF(b)) {
   cudacomplex result = {{ a.value.x / b, a.value.y / b }};
   return result;
}

// divide scalar by complex
HOSTDEVICE cudacomplex operator/(const flt_doub REF(a), const cudacomplex REF(b)) {
   flt_doub tmp = ( b.value.x * b.value.x + b.value.y * b.value.y );
   cudacomplex result = {{ ( a * b.value.x ) / tmp, ( -a * b.value.y ) / tmp }};
   return result;
}

// complex conjugate
HOSTDEVICE cudacomplex operator~(const cudacomplex REF(a)) {
   cudacomplex result = {{ a.value.x, -a.value.y }};
   return result;
}

// a possible alternative to a cudacomplex constructor
HOSTDEVICE cudacomplex make_cudacomplex(flt_doub a, flt_doub b)
{
    cudacomplex res;
    res.real() = a;
    res.imag() = b;
    return res;
}
/*
namespace constants
{
    const _cudacomplex zero = make_cudacomplex(0.0f, 0.0f);
    const _cudacomplex one  = make_cudacomplex(1.0f, 0.0f);
    const _cudacomplex I    = make_cudacomplex(0.0f, 1.0f);
};
*/
