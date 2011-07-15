#include <cstring>
#include <cmath>
#include "edge_math.h"

namespace edge {
namespace math {

// Vec2
Vec2::Vec2(double x, double y) {    
    vec[0] = x;
    vec[1] = y;
}

Vec2::Vec2(const double *arr2) {    
    vec[0] = arr2[0];
    vec[1] = arr2[1];
}

Vec2::Vec2(const Vec2 &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
}

Vec2::Vec2(const Vec4 &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
}

Vec2::Vec2(const Vec3 &v) {
    vec[0] = v.vec[0];
    vec[1] = v.vec[1];
}

Vec2 operator+(const Vec2 &op1, const Vec2 &op2) {
    Vec2 result;
    result.vec[0] = op1.vec[0] + op2.vec[0];
    result.vec[1] = op1.vec[1] + op2.vec[1];

    return result;
}

Vec2 operator-(const Vec2 &op1, const Vec2 &op2) {
    Vec2 result;
    result.vec[0] = op1.vec[0] - op2.vec[0];
    result.vec[1] = op1.vec[1] - op2.vec[1];

    return result;
}

double operator*(const Vec2 &op1, const Vec2 &op2) {
    return op1.vec[0]*op2.vec[0] + 
           op1.vec[1]*op2.vec[1];
}

Vec2 operator*(const Vec2 &op1, const double op2) {
    Vec2 result;
    result.vec[0] = op1.vec[0]*op2;
    result.vec[1] = op1.vec[1]*op2;

    return result;
}

Vec2 operator*(const double op1, const Vec2 &op2) {
    Vec2 result;
    result.vec[0] = op2.vec[0]*op1;
    result.vec[1] = op2.vec[1]*op1;

    return result;
}

double Vec2::magnitude() const {
    return sqrt(this->vec[0]*this->vec[0] + 
                this->vec[1]*this->vec[1]);
}

} // namespace math
} // namespace edge