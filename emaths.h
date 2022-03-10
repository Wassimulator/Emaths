/*
    Emaths - Essential Maths library

    includes most 2D geometry related math functions
    Work in progress, more functions will be added as needed

    by Wassimulator
*/

#pragma once
#include <stdio.h>
#include <iostream>

float clamp(float input, float min, float max)
{
    if (input < min)
        return min;
    else if (input > max)
        return max;
    else
        return input;
}
int rand_range(int min, int max)
{
    return min + rand() / (RAND_MAX / (max - min + 1) + 1);
}
inline double random_double()
{
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_double();
}
struct v2
{
    union
    {
        float e[2];
        struct
        {
            float x, y;
        };
    };

    v2() : e{0, 0} {}
    v2(float e0, float e1) : e{e0, e1} {}

    v2 operator-() const { return v2(-e[0], -e[1]); }
    float operator[](int i) const { return e[i]; }
    float &operator[](int i) { return e[i]; }

    v2 &operator+=(const v2 &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        return *this;
    }
    v2 &operator-=(const v2 &v)
    {
        e[0] -= v.e[0];
        e[1] -= v.e[1];
        return *this;
    }

    v2 &operator*=(const float t)
    {
        e[0] *= t;
        e[1] *= t;
        return *this;
    }

    v2 &operator/=(const float t)
    {
        return *this *= 1 / t;
    }

    float length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1];
    }

    float length() const
    {
        return sqrt(length_squared());
    }

    v2 normalize()
    {
        float lengf = length();

        if (lengf > 0)
            return v2(x / lengf, y / lengf);
        else
            return v2(0, 0);
    }

    float dot(v2 V)
    {
        float result = x * V.x + y * V.y;
        return result;
    }
    v2 perpendicular()
    {
        return v2(y, -x);
    }
    float perpdot(v2 V)
    {
        float result = x * V.y - y * V.x;
        return result;
    }
    float cross(v2 V)
    {
        float result = x * V.y - y * V.x;
        return result;
    }

    v2 hadamard(v2 V)
    {
        v2 result(x * V.x, y * V.y);
        return result;
    }

    float angle(v2 V) // returns signed angle in radians
    {
        return atan2(x * V.y - y * V.x, x * V.x + y * V.y);
    }
};
v2 operator-(v2 a, v2 b)
{
    v2 tojesus(a.x - b.x, a.y - b.y);
    return tojesus;
}

v2 operator+(v2 a, v2 b)
{
    v2 tojesus(a.x + b.x, a.y + b.y);
    return tojesus;
}
v2 operator-(v2 a, float b)
{
    v2 tojesus(a.x - b, a.y - b);
    return tojesus;
}

v2 operator+(v2 a, float b)
{
    v2 tojesus(a.x + b, a.y + b);
    return tojesus;
}

v2 operator*(v2 a, float b)
{
    v2 tojesus(a.x * b, a.y * b);
    return tojesus;
}
v2 operator*(v2 a, v2 b)
{
    v2 tojesus(a.x * b.x, a.y * b.y);
    return tojesus;
}

v2 operator*(float b, v2 a)
{
    v2 tojesus(a.x * b, a.y * b);
    return tojesus;
}

v2 operator/(v2 a, float b)
{
    v2 tojesus(a.x / b, a.y / b);
    return tojesus;
}
v2 operator/(v2 a, v2 b)
{
    v2 tojesus(a.x / b.x, a.y / b.y);
    return tojesus;
}

float v2_distance_2Points(v2 A, v2 B)
{
    return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

v2 unitvec_AtoB(v2 A, v2 B)
{
    float n = v2_distance_2Points(A, B);
    return ((B - A) / n);
}

float signed_angle_v2(v2 A, v2 B)
{
    return atan2(A.x * B.y - A.y * B.x, A.x * B.x + A.y * B.y);
}

struct v3
{
    union
    {
        float e[3];
        struct
        {
            float x, y, z;
        };
    };

    v3() : e{0, 0, 0} {}
    v3(float e0, float e1, float e2) : e{e0, e1, e2} {}

    v3 operator-() const { return v3(-e[0], -e[1], -e[2]); }
    float operator[](int i) const { return e[i]; }
    float &operator[](int i) { return e[i]; }

    v3 &operator+=(const v3 &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    v3 &operator*=(const float t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    v3 &operator/=(const float t)
    {
        return *this *= 1 / t;
    }

    float length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    float length() const
    {
        return sqrt(length_squared());
    }

    v3 normalize()
    {
        float lengf = length();
        v3 thomas(x / lengf, y / lengf, z / lengf);
        return thomas;
    }

    float dot(v3 V)
    {
        float result = x * V.x + y * V.y + z * V.z;
        return result;
    }

    v3 hadamard(v3 V)
    {
        v3 result(x * V.x, y * V.y, z * V.z);
        return result;
    }
};

v3 operator-(v3 a, v3 b)
{
    v3 tojesus(a.x - b.x, a.y - b.y, a.z - b.z);
    return tojesus;
}

v3 operator+(v3 a, v3 b)
{
    v3 tojesus(a.x + b.x, a.y + b.y, a.z + b.z);
    return tojesus;
}

v3 operator*(v3 a, float b)
{
    v3 tojesus(a.x * b, a.y * b, a.z * b);
    return tojesus;
}

v3 operator*(float b, v3 a)
{
    v3 tojesus(a.x * b, a.y * b, a.z * b);
    return tojesus;
}

v3 operator/(v3 a, float b)
{
    v3 tojesus(a.x / b, a.y / b, a.z / b);
    return tojesus;
}

// Type aliases for vec3
using point3 = v3; // 3D point
// using color = v3;  // RGB color

v2 Rotate2D(v2 P, float sine, float cosine)
{
    return v2(v2(cosine, -sine).dot(P), v2(sine, cosine).dot(P));
}
v2 Rotate2D(v2 P, float Angle)
{
    float sine = sin(Angle);
    float cosine = cos(Angle);
    return Rotate2D(P, sine, cosine);
}
v2 Rotate2D(v2 p, v2 axis, float angle)
{
    // POINT rotate_point(float cx,float cy,float angle,POINT p)
    float s = sin(angle);
    float c = cos(angle);

    // translate point back to origin:
    p.x -= axis.x;
    p.y -= axis.y;

    // rotate point
    float xnew = p.x * c - p.y * s;
    float ynew = p.x * s + p.y * c;

    // translate point back:
    p.x = xnew + axis.x;
    p.y = ynew + axis.y;
    return p;
}

v2 Reflection2D(v2 P, v2 N)
{
    return P - 2 * N.dot(P) * N;
}

bool PointInRectangle(v2 P, v2 A, v2 B, v2 C)
{
    v2 M = P;
    v2 AB = B - A;
    v2 BC = C - B;
    v2 AM = M - A;
    v2 BM = M - B;
    if (0 <= AB.dot(AM) && AB.dot(AM) <= AB.dot(AB) &&
        0 <= BC.dot(BM) && BC.dot(BM) <= BC.dot(BC))
        return true;
    else
        return false;
}

int sign(float x)
{
    if (x > 0)
        return 1;
    else if (x < 0)
        return -1;
    else
        return 0;
}

static float dot(v2 A, v2 B)
{
    return A.x * B.x + A.y * B.y;
}

static float perpdot(v2 A, v2 B)
{
    return A.x * B.y - A.y * B.x;
}

static bool operator==(v2 A, v2 B)
{
    return A.x == B.x && A.y == B.y;
}

// POINT rotate_point(v2 around, float angle, POINT p)
// {
//     float cx = around.x;
//     float cy = around.y;
//     float s = sin(angle);
//     float c = cos(angle);

//     // translate point back to origin:
//     p.x -= cx;
//     p.y -= cy;

//     // rotate point
//     float xnew = p.x * c - p.y * s;
//     float ynew = p.x * s + p.y * c;

//     // translate point back:
//     p.x = xnew + cx;
//     p.y = ynew + cy;
//     return p;
// }