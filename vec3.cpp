#include <cmath>
#include "edge_math.h"

namespace edge {
namespace math {

Vec3::Vec3(double x, double y, double z) {
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
}

Vec3::Vec3(const double *arr3) {
    vec[0] = arr3[0];
    vec[1] = arr3[1];
    vec[2] = arr3[2];
    vec[3] = arr3[3];
}

Vec3::Vec3(const Vec3 &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
    vec[2] = v.vec[2];
}


Vec3 operator+(const Vec3 &op1, const Vec3 &op2) {
    Vec3 result;
    result.vec[0] = op1.vec[0] + op2.vec[0];
    result.vec[1] = op1.vec[1] + op2.vec[1];
    result.vec[2] = op1.vec[2] + op2.vec[2];

    return result;
}

Vec3 operator-(const Vec3 &op1, const Vec3 &op2) {
    Vec3 result;
    result.vec[0] = op1.vec[0] - op2.vec[0];
    result.vec[1] = op1.vec[1] - op2.vec[1];
    result.vec[2] = op1.vec[2] - op2.vec[2];

    return result;
}

double operator*(const Vec3 &op1, const Vec3 &op2) {
    return op1.vec[0]*op2.vec[0] + 
           op1.vec[1]*op2.vec[1] + 
           op1.vec[2]*op2.vec[2];
}

Vec3 operator*(const Vec3 &op1, const double op2) {
    Vec3 result;
    result.vec[0] = op1.vec[0]*op2;
    result.vec[1] = op1.vec[1]*op2;
    result.vec[2] = op1.vec[2]*op2;

    return result;
}

Vec3 operator*(const double op1, const Vec3 &op2) {
    Vec3 result;
    result.vec[0] = op2.vec[0]*op1;
    result.vec[1] = op2.vec[1]*op1;
    result.vec[2] = op2.vec[2]*op1;

    return result;
}


double Vec3::magnitude() const {
    return sqrt(this->vec[0]*this->vec[0] + 
                this->vec[1]*this->vec[1] + 
                this->vec[2]*this->vec[2]);
}

Vec3 Vec3::crossProduct(const Vec3 &op) {
    Vec3 result;
    result.vec[0] = this->vec[1]*op.vec[2] - this->vec[2]*op.vec[1];
    result.vec[1] = -this->vec[0]*op.vec[2] + this->vec[2]*op.vec[0];
    result.vec[2] = this->vec[0]*op.vec[1] - this->vec[1]*op.vec[0];
    return result;
}

} // namespace math
} // namespace edge