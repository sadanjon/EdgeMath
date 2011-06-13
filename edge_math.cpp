#include <cmath>
#include "edge_math.h"

using namespace edge::math;

// Vec4
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

// Vec3
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

// Vec2
Vec2::Vec2(double x, double y) {
    vec[0] = x;
    vec[1] = y;
}

Vec2::Vec2(const double *arr2) {
    vec[0] = arr2[0];
    vec[1] = arr2[1];
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

// Mat44
Mat44::Mat44(const double *arr16) {
    for (int i = 15; i >= 0; --i)
        mat[i] = arr16[i];
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

#define M(i,j) op1.mat[4*i + j]
Vec4 operator*(const Mat44 &op1, const Vec4 &op2) {
    Vec4 result;
    
    result.vec[0] = M(0,0)*op2.vec[0] +
                    M(0,1)*op2.vec[1] +
                    M(0,2)*op2.vec[2] +
                    M(0,3)*op2.vec[3];

    result.vec[1] = M(1,0)*op2.vec[0] +
                    M(1,1)*op2.vec[1] +
                    M(1,2)*op2.vec[2] +
                    M(1,3)*op2.vec[3];

    result.vec[2] = M(2,0)*op2.vec[0] +
                    M(2,1)*op2.vec[1] +
                    M(2,2)*op2.vec[2] +
                    M(2,3)*op2.vec[3];

    result.vec[3] = M(3,0)*op2.vec[0] +
                    M(3,1)*op2.vec[1] +
                    M(3,2)*op2.vec[2] +
                    M(3,3)*op2.vec[3];

    return result;
}
#undef M

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

void Mat44::insertCol(int col, const Vec4 &vec) {
    insertCol(col, vec.vec);
}

void Mat44::insertCol(int col, const double *vec) {
    col *= 4;
    mat[col] = vec[0];
    mat[col + 1] = vec[1];
    mat[col + 2] = vec[2];
    mat[col + 3] = vec[3];
}

#define M(i,j) mat[i*j + j]
void Mat44::transpose() {
    for (int i = 3; i >= 0; --i) {
        double temp;

        temp = M(i, 0);
        M(i, 0) = M(0, i);
        M(0, i) = temp;

        temp = M(i, 1);
        M(i, 1) = M(1, i);
        M(1, i) = temp;

        temp = M(i, 2);
        M(i, 2) = M(2, i);
        M(2, i) = temp;

        temp = M(i, 3);
        M(i, 3) = M(3, i);
        M(3, i) = temp;
    }
}


Mat33 Mat44::getMinor(int i, int j) const {
    Mat33 result;

    int res_r = 0, res_c = 0;
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            if (r != i && c != j) {
                result.mat[res_r*3 + res_c] = this->mat[r*4 + c];
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

void Mat44::inverse() {
    Mat44 inverseMat = getInverse();
    for (int i = 15; i >= 0; --i)
        mat[i] = inverseMat.mat[i];

}
#undef M

// Mat33
Mat33::Mat33() {
}

Mat33::Mat33(DOUBLE *arr9) {
    for (int i = 8; i >= 0; --i)
        mat[i] = arr9[i];
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

#define M(i,j) op1.mat[3*i + j]
Vec3 operator*(const Mat33 &op1, const Vec3 &op2) {
    Vec3 result;
    
    result.vec[0] = M(0,0)*op2.vec[0] +
                    M(0,1)*op2.vec[1] +
                    M(0,2)*op2.vec[2];

    result.vec[1] = M(1,0)*op2.vec[0] +
                    M(1,1)*op2.vec[1] +
                    M(1,2)*op2.vec[2];

    result.vec[2] = M(2,0)*op2.vec[0] +
                    M(2,1)*op2.vec[1] +
                    M(2,2)*op2.vec[2];

    return result;
}
#undef M

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

void Mat33::insertCol(int col, const Vec3 &vec) {
    insertCol(col, vec.vec);
}

void Mat33::insertCol(int col, const double *vec) {
    col *= 3;
    mat[col] = vec[0];
    mat[col + 1] = vec[1];
    mat[col + 2] = vec[2];
}

#define M(i,j) mat[i*j + j]
void Mat33::transpose() {
    for (int i = 2; i >= 0; --i) {
        double temp;

        temp = M(i, 0);
        M(i, 0) = M(0, i);
        M(0, i) = temp;

        temp = M(i, 1);
        M(i, 1) = M(1, i);
        M(1, i) = temp;

        temp = M(i, 2);
        M(i, 2) = M(2, i);
        M(2, i) = temp;
        
    }
}
#undef M

Mat22 Mat33::getMinor(int i, int j) const{
    Mat22 result;

    int res_r = 0, res_c = 0;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            if (r != i && c != j) {
                result.mat[res_r*2 + res_c] = this->mat[r*3 + c];
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

#define M(i,j) mat[i*4 + j]
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

void Mat33::inverse() {
    Mat33 inverseMat = getInverse();
    for (int i = 8; i >= 0; --i)
        mat[i] = inverseMat[i];
}

#undef M

// Mat22
Mat22::Mat22() {
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

#define M(i,j) op1.mat[2*i + j]
Vec2 operator*(const Mat22 &op1, const Vec2 &op2) {
    Vec2 result;
    
    result.vec[0] = M(0,0)*op2.vec[0] +
                    M(0,1)*op2.vec[1];

    result.vec[1] = M(1,0)*op2.vec[0] +
                    M(1,1)*op2.vec[1];    

    return result;
}
#undef M

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

void Mat22::insertCol(int col, const Vec2 &vec) {
    insertCol(col, vec.vec);
}

void Mat22::insertCol(int col, const double *vec) {
    col *= 2;
    mat[col] = vec[0];
    mat[col + 1] = vec[1];    
}


void Mat22::transpose() {
    double temp;

    temp = mat[2];
    mat[1] = mat[2];
    mat[2] = temp;
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

void Mat22::inverse() {
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
}

