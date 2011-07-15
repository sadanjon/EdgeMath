#include <cstring>
#include "edge_math.h"

namespace edge {
namespace math {

#define M(i,j) mat[j*4 + i]

Mat44::Mat44(const double *arr16) {
    memcpy(this->mat, arr16, 16*sizeof(double));    
}

Mat44::Mat44(const Mat44 &m) {
    memcpy(this->mat, m.mat, 16*sizeof(double));    
}

Mat44 operator+(const Mat44 &op1, const Mat44 &op2) {
    Mat44 result;
    for (int i = 15; i >= 0; --i) 
        result.mat[i] = op1.mat[i] + op2.mat[i];
    
    return result;
}
Mat44 operator-(const Mat44 &op1, const Mat44 &op2) {
    Mat44 result;
    for (int i = 15; i >= 0; --i) 
        result.mat[i] = op1.mat[i] - op2.mat[i];
    
    return result;
}
Mat44 operator*(const Mat44 &op1, const Mat44 &op2) {
    Mat44 result;
    for (int i = 3; i >= 0; --i) {
        int col = 4*i;
        result.mat[col] = op1.mat[0]*op2.mat[col] +
                            op1.mat[4]*op2.mat[col + 1] +
                            op1.mat[8]*op2.mat[col + 2] +
                            op1.mat[12]*op2.mat[col + 3];

        result.mat[col + 1] = op1.mat[1]*op2.mat[col] +
                            op1.mat[5]*op2.mat[col + 1] +
                            op1.mat[9]*op2.mat[col + 2] +
                            op1.mat[13]*op2.mat[col + 3];

        result.mat[col + 2] = op1.mat[2]*op2.mat[col] +
                            op1.mat[6]*op2.mat[col + 1] +
                            op1.mat[10]*op2.mat[col + 2] +
                            op1.mat[14]*op2.mat[col + 3];

        result.mat[col + 3] = op1.mat[3]*op2.mat[col] +
                            op1.mat[7]*op2.mat[col + 1] +
                            op1.mat[11]*op2.mat[col + 2] +
                            op1.mat[15]*op2.mat[col + 3];
    }

    return result;
}


Vec4 operator*(const Mat44 &op1, const Vec4 &op2) {
    Vec4 result;
    
    result.vec[0] = op1.M(0,0)*op2.vec[0] +
                    op1.M(0,1)*op2.vec[1] +
                    op1.M(0,2)*op2.vec[2] +
                    op1.M(0,3)*op2.vec[3];

    result.vec[1] = op1.M(1,0)*op2.vec[0] +
                    op1.M(1,1)*op2.vec[1] +
                    op1.M(1,2)*op2.vec[2] +
                    op1.M(1,3)*op2.vec[3];

    result.vec[2] = op1.M(2,0)*op2.vec[0] +
                    op1.M(2,1)*op2.vec[1] +
                    op1.M(2,2)*op2.vec[2] +
                    op1.M(2,3)*op2.vec[3];

    result.vec[3] = op1.M(3,0)*op2.vec[0] +
                    op1.M(3,1)*op2.vec[1] +
                    op1.M(3,2)*op2.vec[2] +
                    op1.M(3,3)*op2.vec[3];

    return result;
}


Mat44 operator*(const Mat44 &op1, const double op2) {
    Mat44 result;
    for (int i = 15; i >= 0; --i)
        result.mat[i] = op1.mat[i] * op2;

    return result;    
}

Mat44 operator*(const double op1, const Mat44 &op2) {
    Mat44 result;
    for (int i = 15; i >= 0; --i)
        result.mat[i] = op2.mat[i] * op1;

    return result;    
}

Mat44 &Mat44::insertCol(int col, const Vec4 &vec) {
    insertCol(col, vec.vec);

    return *this;
}

Mat44 &Mat44::insertCol(int col, const double *vec) {
    col *= 4;
    mat[col] = vec[0];
    mat[col + 1] = vec[1];
    mat[col + 2] = vec[2];
    mat[col + 3] = vec[3];

    return *this;
}

Mat44 Mat44::getTranspose() const {
    Mat44 m(*this);
    m.transpose();
    return m;
}


Mat44 &Mat44::transpose() {
    double temp;

    temp = M(1, 0);
    M(1, 0) = M(0, 1);
    M(0, 1) = temp;

    temp = M(2, 0);
    M(2, 0) = M(0, 2);
    M(0, 2) = temp;

    temp = M(3, 0);
    M(3, 0) = M(0, 3);
    M(0, 3) = temp;

    temp = M(3, 1);
    M(3, 1) = M(1, 3);
    M(1, 3) = temp;

    temp = M(3, 2);
    M(3, 2) = M(2, 3);
    M(2, 3) = temp;

    temp = M(2, 1);
    M(2, 1) = M(1, 2);
    M(1, 2) = temp;


    return *this;
}


Mat33 Mat44::getMinor(int i, int j) const {
    Mat33 result;

    int res_r = 0, res_c = 0;
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            if (r != i && c != j) {
                result.mat[res_c*3 + res_r] = this->M(r, c);
                ++res_c;
                if (res_c == 3) {
                    res_c = 0;
                    ++res_r;
                }
            }
        }
    }

    return result;
}

double Mat44::getDeterminant() const {  
    Mat33 m00 = getMinor(0,0);
    Mat33 m01 = getMinor(0,1);
    Mat33 m02 = getMinor(0,2);
    Mat33 m03 = getMinor(0,3);

    return M(0,0)*m00.getDeterminant() - M(0,1)*m01.getDeterminant() +
        M(0,2)*m02.getDeterminant() - M(0,3)*m03.getDeterminant();
}

Mat44 Mat44::getInverse() const {
    Mat44 result;
    double det_recp = 1/getDeterminant();

    for (int i = 3; i >= 0; --i) {        
        double minor_det = getMinor(0,i).getDeterminant();        
        result.M(i,0) = det_recp*minor_det*((i+0)%2 == 0 ? 1: -1);

        minor_det = getMinor(1,i).getDeterminant();
        result.M(i,1) = det_recp*minor_det*((i+1)%2 == 0 ? 1: -1);

        minor_det = getMinor(2,i).getDeterminant();
        result.M(i,2) = det_recp*minor_det*((i+2)%2 == 0 ? 1: -1);

        minor_det = getMinor(3,i).getDeterminant();
        result.M(i,3) = det_recp*minor_det*((i+3)%2 == 0 ? 1: -1);
    }

    return result;

}

Mat44 &Mat44::inverse() {
    Mat44 inverseMat = getInverse();
    memcpy(mat, inverseMat.mat, 16*sizeof(double));    

    return *this;
}


} // namespace math
} // namespace edge