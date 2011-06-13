#include <cmath>
#include "edge_math.h"


namespace edge {
namespace math {

Vec4::Vec4(double x, double y, double z, double w) {
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
    vec[3] = w;
}

Vec4::Vec4(const double *arr4) {
    vec[0] = arr4[0];
    vec[1] = arr4[1];
    vec[2] = arr4[2];
    vec[3] = arr4[3];
}

Vec4::Vec4(const Vec4 &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
    vec[2] = v.vec[2];
    vec[3] = v.vec[3];
}


Vec4 operator+(const Vec4 &op1, const Vec4 &op2) {
    Vec4 result;
    result.vec[0] = op1.vec[0] + op2.vec[0];
    result.vec[1] = op1.vec[1] + op2.vec[1];
    result.vec[2] = op1.vec[2] + op2.vec[2];
    result.vec[3] = op1.vec[3] + op2.vec[3];

    return result;
}

Vec4 operator-(const Vec4 &op1, const Vec4 &op2) {
    Vec4 result;
    result.vec[0] = op1.vec[0] - op2.vec[0];
    result.vec[1] = op1.vec[1] - op2.vec[1];
    result.vec[2] = op1.vec[2] - op2.vec[2];
    result.vec[3] = op1.vec[3] - op2.vec[3];

    return result;
}

double operator*(const Vec4 &op1, const Vec4 &op2) {
    return op1.vec[0]*op2.vec[0] + 
           op1.vec[1]*op2.vec[1] + 
           op1.vec[2]*op2.vec[2] + 
           op1.vec[3]*op2.vec[3];
}

Vec4 operator*(const Vec4 &op1, const double op2) {
    Vec4 result;
    result.vec[0] = op1.vec[0]*op2;
    result.vec[1] = op1.vec[1]*op2;
    result.vec[2] = op1.vec[2]*op2;
    result.vec[3] = op1.vec[3]*op2;

    return result;
}

Vec4 operator*(const double op1, const Vec4 &op2) {
    Vec4 result;
    result.vec[0] = op2.vec[0]*op1;
    result.vec[1] = op2.vec[1]*op1;
    result.vec[2] = op2.vec[2]*op1;
    result.vec[3] = op2.vec[3]*op1;

    return result;
}

double Vec4::magnitude() const {
    return sqrt(this->vec[0]*this->vec[0] + 
                this->vec[1]*this->vec[1] + 
                this->vec[2]*this->vec[2] + 
                this->vec[3]*this->vec[3]);
}

} // namespace math
} // namespace edge