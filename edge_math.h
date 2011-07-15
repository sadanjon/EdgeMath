#pragma once

namespace edge {
namespace math {
    // forward declarations
    class Mat44;
    class Mat33;
    class Mat22;
    class Vec4;
    class Vec3;
    class Vec2;   

    class Vec4 {
    public:
        Vec4(double x = 0, double y = 0, double z = 0, double w = 1);  
        Vec4(const double *arr4);
        Vec4(const Vec4 &v);
        Vec4(const Vec3 &v);
        Vec4(const Vec2 &v);
        double vec[4];

        double magnitude() const;

        friend Vec4 operator+(const Vec4 &op1, const Vec4 &op2);
        friend Vec4 operator-(const Vec4 &op1, const Vec4 &op2);
        friend double operator*(const Vec4 &op1, const Vec4 &op2);
        friend Vec4 operator*(const Vec4 &op1, const double op2);
        friend Vec4 operator*(const double op1, const Vec4 &op2);
    };    

    class Vec3 {
    public:
        Vec3(double x = 0, double y = 0, double z = 0);
        Vec3(const double *arr3);
        Vec3(const Vec3 &v);
        Vec3(const Vec4 &v);
        Vec3(const Vec2 &v);
        double vec[3];

        double magnitude() const;
        static Vec3 crossProduct(const Vec3 &op1, const Vec3 &op2);

        friend Vec3 operator+(const Vec3 &op1, const Vec3 &op2);
        friend Vec3 operator-(const Vec3 &op1, const Vec3 &op2);
        friend double operator*(const Vec3 &op1, const Vec3 &op2); 
        friend Vec3 operator*(const Vec3 &op1, const double op2);
        friend Vec3 operator*(const double op1, const Vec3 &op2);
 
    };

    class Vec2 {
    public:
        Vec2(double x = 0, double y = 0);
        Vec2(const double *arr2);
        Vec2(const Vec2 &v);
        Vec2(const Vec4 &v);
        Vec2(const Vec3 &v);
        double vec[2];

        double magnitude() const;

        friend Vec2 operator+(const Vec2 &op1, const Vec2 &op2);
        friend Vec2 operator-(const Vec2 &op1, const Vec2 &op2);
        friend double operator*(const Vec2 &op1, const Vec2 &op2);
        friend Vec2 operator*(const Vec2 &op1, const double op2);
        friend Vec2 operator*(const double op1, const Vec2 &op2);
    };

    

    class Mat44 {
    public:
        // initializes with junk
        Mat44() {};
        Mat44(const double *arr16);   
        Mat44(const Mat44 &m);

        double mat[16];

        // col is between 0 and 3
        Mat44 &insertCol(int col, const Vec4 &vec);
        Mat44 &insertCol(int col, const double *vec4);
        Mat44 &transpose();
        // assumes det != 0. will divide by zero else!
        Mat44 &inverse();
        
        double getDeterminant() const;        
        Mat44 getTranspose() const;
        // assumes det != 0. will divide by zero else!
        Mat44 getInverse() const;
        Mat33 getMinor(int i, int j) const;

        static Mat44 getIdentity();
        

        friend Mat44 operator+(const Mat44 &op1, const Mat44 &op2);
        friend Mat44 operator-(const Mat44 &op1, const Mat44 &op2);
        friend Mat44 operator*(const Mat44 &op1, const Mat44 &op2);
        friend Vec4 operator*(const Mat44 &op1, const Vec4 &op2);
        friend Mat44 operator*(const Mat44 &op1, const double op2);
        friend Mat44 operator*(const double op1, const Mat44 &op2);

    };

    class Mat33 {
     public:
        Mat33() {}
        Mat33(const double *arr9);
        Mat33(const Mat33 &m);
        Mat33(const Mat33 *m);
        double mat[9];

        // col is between 0 and 2
        Mat33 &insertCol(int col, const Vec3 &vec);
        Mat33 &insertCol(int col, const double *vec3);
        Mat33 &transpose();
        // assumes det != 0. will divide by zero else!
        Mat33 &inverse();
        
        double getDeterminant() const;        
        Mat33 getTranspose() const;
        // assumes det != 0. will divide by zero else!
        Mat33 getInverse() const;
        Mat22 getMinor(int i, int j) const;

        static Mat33 getIdentity();

        friend Mat33 operator+(const Mat33 &op1, const Mat33 &op2);
        friend Mat33 operator-(const Mat33 &op1, const Mat33 &op2);
        friend Mat33 operator*(const Mat33 &op1, const Mat33 &op2);
        friend Vec3 operator*(const Mat33 &op1, const Vec3 &op2);
        friend Mat33 operator*(const Mat33 &op1, const double op2);
        friend Mat33 operator*(const double op1, const Mat33 &op2);
    };

    class Mat22 {
     public:
         Mat22() {}
        Mat22(double *arr4);
        Mat22(const Mat22 &m);
        double mat[4];

         // col is between 0 and 1
        Mat22 &insertCol(int col, const Vec2 &vec);
        Mat22 &insertCol(int col, const double *vec2);
        Mat22 &transpose();
        // assumes det != 0. will divide by zero else!
        Mat22 &inverse();
        
        double getDeterminant() const;        
        Mat22 getTranspose() const;
        // assumes det != 0. will divide by zero else!
        Mat22 getInverse() const;        

        static Mat22 getIdentity();

        friend Mat22 operator+(const Mat22 &op1, const Mat22 &op2);
        friend Mat22 operator-(const Mat22 &op1, const Mat22 &op2);
        friend Mat22 operator*(const Mat22 &op1, const Mat22 &op2);
        friend Vec2 operator*(const Mat22 &op1, const Vec2 &op2);
        friend Mat22 operator*(const Mat22 &op1, const double op2);
        friend Mat22 operator*(const double op1, const Mat22 &op2);
    };

    int *nextPerm4(bool reset);
    
        
} // namespace math
} // namespace edge