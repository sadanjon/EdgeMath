#include "edge_math.h"

namespace edge {
namespace math {

#define M(i,j) mat[2*j + i]

Mat22::Mat22(double *arr4) {
    mat[0] = arr4[0];
    mat[1] = arr4[1];
    mat[2] = arr4[2];
    mat[3] = arr4[3];
}

Mat22::Mat22(const Mat22 &m) {
    mat[0] = m.mat[0];
    mat[1] = m.mat[1];
    mat[2] = m.mat[2];
    mat[3] = m.mat[3];
}


Mat22 operator+(const Mat22 &op1, const Mat22 &op2) {
    Mat22 result;
    for (int i = 3; i >= 0; --i) 
        result.mat[i] = op1.mat[i] + op2.mat[i];
    
    return result;
}
Mat22 operator-(const Mat22 &op1, const Mat22 &op2) {
    Mat22 result;
    for (int i = 3; i >= 0; --i) 
        result.mat[i] = op1.mat[i] - op2.mat[i];
    
    return result;
}

Mat22 operator*(const Mat22 &op1, const Mat22 &op2) {
    Mat22 result;

    result.mat[0] = op1.mat[0]*op2.mat[0] +
        op1.mat[2]*op2.mat[1];


    result.mat[1] = op1.mat[1]*op2.mat[0] +
        op1.mat[3]*op2.mat[1];

    result.mat[2] = op1.mat[0]*op2.mat[2] +
        op1.mat[2]*op2.mat[3];


    result.mat[3] = op1.mat[1]*op2.mat[2] +
        op1.mat[3]*op2.mat[3];

    return result;
}


Vec2 operator*(const Mat22 &op1, const Vec2 &op2) {
    Vec2 result;
    
    result.vec[0] = op1.M(0,0)*op2.vec[0] +
                    op1.M(0,1)*op2.vec[1];

    result.vec[1] = op1.M(1,0)*op2.vec[0] +
                    op1.M(1,1)*op2.vec[1];    

    return result;
}


Mat22 operator*(const Mat22 &op1, const double op2) {
    Mat22 result;
    for (int i = 3; i >= 0; --i)
        result.mat[i] = op1.mat[i] * op2;

    return result;    
}

Mat22 operator*(const double op1, const Mat22 &op2) {
    Mat22 result;
    for (int i = 3; i >= 0; --i)
        result.mat[i] = op2.mat[i] * op1;

    return result;    
}

Mat22 &Mat22::insertCol(int col, const Vec2 &vec) {
    insertCol(col, vec.vec);

    return *this;
}

Mat22 &Mat22::insertCol(int col, const double *vec) {
    col *= 2;
    mat[col] = vec[0];
    mat[col + 1] = vec[1];    

    return *this;
}

Mat22 Mat22::getTranspose() const {
    return Mat22(*this).transpose();
}

Mat22 &Mat22::transpose() {
    double temp;

    temp = mat[2];
    mat[1] = mat[2];
    mat[2] = temp;

    return *this;
}

double Mat22::getDeterminant() const {
    return mat[0]*mat[3] - mat[2]*mat[1];
}

Mat22 Mat22::getInverse() const {
    double arr[4];
    double det_recp = 1/getDeterminant();

    arr[0] = det_recp*mat[3];
    arr[1] = -det_recp*mat[1];
    arr[2] = -det_recp*mat[2];
    arr[3] = det_recp*mat[0];

    return Mat22(arr);
}

Mat22 &Mat22::inverse() {
    double arr[4];
    double det_recp = 1/getDeterminant();

    arr[0] = det_recp*mat[3];
    arr[1] = -det_recp*mat[1];
    arr[2] = -det_recp*mat[2];
    arr[3] = det_recp*mat[0];

    mat[0] = arr[0];
    mat[1] = arr[1];
    mat[2] = arr[2];
    mat[3] = arr[3];

    return *this;
}

} // namespace math
} // namespace edge