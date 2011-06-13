#include <cstring>
#include "edge_math.h"

namespace edge {
namespace math {

#define M(i,j) mat[j*3 + i]

Mat33::Mat33(const double *arr9) {
    memcpy(mat, arr9, 9*sizeof(double));    
}

Mat33::Mat33(const Mat33 &m) {
    memcpy(mat, m.mat, 9*sizeof(double));        
}

Mat33 operator+(const Mat33 &op1, const Mat33 &op2) {
    Mat33 result;
    for (int i = 8; i >= 0; --i) 
        result.mat[i] = op1.mat[i] + op2.mat[i];
    
    return result;
}
Mat33 operator-(const Mat33 &op1, const Mat33 &op2) {
    Mat33 result;
    for (int i = 8; i >= 0; --i) 
        result.mat[i] = op1.mat[i] - op2.mat[i];
    
    return result;
}
Mat33 operator*(const Mat33 &op1, const Mat33 &op2) {
    Mat33 result;
    for (int i = 2; i >= 0; --i) {
        int col = 3*i;
        result.mat[col] = op1.mat[0]*op2.mat[col] +
                            op1.mat[3]*op2.mat[col + 1] +
                            op1.mat[6]*op2.mat[col + 2];

        result.mat[col + 1] = op1.mat[1]*op2.mat[col] +
                            op1.mat[4]*op2.mat[col + 1] +
                            op1.mat[7]*op2.mat[col + 2];                            

        result.mat[col + 2] = op1.mat[2]*op2.mat[col] +
                            op1.mat[5]*op2.mat[col + 1] +
                            op1.mat[8]*op2.mat[col + 2];                        
    }

    return result;
}

Vec3 operator*(const Mat33 &op1, const Vec3 &op2) {
    Vec3 result;
    
    result.vec[0] = op1.M(0,0)*op2.vec[0] +
                    op1.M(0,1)*op2.vec[1] +
                    op1.M(0,2)*op2.vec[2];

    result.vec[1] = op1.M(1,0)*op2.vec[0] +
                    op1.M(1,1)*op2.vec[1] +
                    op1.M(1,2)*op2.vec[2];

    result.vec[2] = op1.M(2,0)*op2.vec[0] +
                    op1.M(2,1)*op2.vec[1] +
                    op1.M(2,2)*op2.vec[2];

    return result;
}

Mat33 operator*(const Mat33 &op1, const double op2) {
    Mat33 result;
    for (int i = 8; i >= 0; --i)
        result.mat[i] = op1.mat[i] * op2;

    return result;    
}

Mat33 operator*(const double op1, const Mat33 &op2) {
    Mat33 result;
    for (int i = 8; i >= 0; --i)
        result.mat[i] = op2.mat[i] * op1;

    return result;    
}

Mat33 &Mat33::insertCol(int col, const Vec3 &vec) {
    insertCol(col, vec.vec);

    return *this;
}

Mat33 &Mat33::insertCol(int col, const double *vec) {
    col *= 3;    
    mat[col] = vec[0];
    mat[col + 1] = vec[1];
    mat[col + 2] = vec[2];

    return *this;
}

Mat33 Mat33::getTranspose() const {        
    Mat33 m(*this);
    m.transpose();
    return m;
}

Mat33 &Mat33::transpose() {    
    double temp;

    temp = M(1, 0);
    M(1, 0) = M(0, 1);
    M(0, 1) = temp;

    temp = M(2, 0);
    M(2, 0) = M(0, 2);
    M(0, 2) = temp;

    temp = M(2, 1);
    M(2, 1) = M(1, 2);
    M(1, 2) = temp;


    return *this;
}

Mat22 Mat33::getMinor(int i, int j) const{
    Mat22 result;

    int res_r = 0, res_c = 0;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            if (r != i && c != j) {
                result.mat[res_r*2 + res_c] = this->M(r,c);                
                ++res_c;
                if (res_c == 2) {
                    res_c = 0;
                    ++res_r;
                }
            }
        }
    }

    return result;
}


double Mat33::getDeterminant() const {
    Mat22 m00 = getMinor(0,0);
    Mat22 m01 = getMinor(0,1);
    Mat22 m02 = getMinor(0,2);

    return M(0,0)*m00.getDeterminant() - M(0,1)*m01.getDeterminant() + M(0,2)*m02.getDeterminant();
}

Mat33 Mat33::getInverse() const {
    Mat33 result;
    double det_recp = 1/getDeterminant();

    for (int i = 2; i >= 0; --i) {        
        double minor_det = getMinor(0,i).getDeterminant();        
        result.M(i,0) = det_recp*minor_det*((i+0)%2 == 0 ? 1: -1);

        minor_det = getMinor(1,i).getDeterminant();
        result.M(i,1) = det_recp*minor_det*((i+1)%2 == 0 ? 1: -1);

        minor_det = getMinor(2,i).getDeterminant();
        result.M(i,2) = det_recp*minor_det*((i+2)%2 == 0 ? 1: -1);
    }

    return result;
}

Mat33 &Mat33::inverse() {
    Mat33 inverseMat = getInverse();
    memcpy(mat, inverseMat.mat, 9*sizeof(double));    

    return *this;
}


} // namespace math
} // namespace edge