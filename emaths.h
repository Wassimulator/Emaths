/*
    Emaths - Essential Maths library

    includes simple 2D geometry related math functions in C/C++
    Work in progress, more functions will be added as needed

    Make your life easier by using namespace Emaths

    v 1.0

    by Wassimulator
*/
#pragma once
#include <stdio.h>

namespace Emaths
{
    float clamp(float input,  float min, float max);
    void  clamp(float *input, float min, float max);
    float clamp(int input,  int min, int max);
    void  clamp(int *input, int min, int max);
    int rand_range(int min, int max);
    int *rand_array(int count);
    uint64_t rand_XOR();
    inline double random_double();
    inline double random_double(double min, double max);
    struct v2;
    struct iv2;
    struct v3;
    struct iv3;
    struct RandState;
    float v2_distance_2Points(v2 A, v2 B);
    v2 unitvec_AtoB(v2 A, v2 B);
    float signed_angle_v2(v2 A, v2 B);
    v2 Rotate2D(v2 P, float sine, float cosine);
    v2 Rotate2D(v2 P, float Angle);
    v2 Rotate2D(v2 p, v2 o, float angle);
    v2 Reflection2D(v2 P, v2 N);
    bool PointInRectangle(v2 P, v2 A, v2 B, v2 C);
    int sign(float x);
    static float dot(v2 A, v2 B);
    static float dot(v3 A, v3 B);
    static float perpdot(v2 A, v2 B);
    static bool operator==(v2 A, v2 B);
    iv2 _iv2(v2 A);
    v2  _v2(iv2 A);
    iv3 _iv3(v3 A);
    v3  _v3(iv3 A);
    uint32_t hash_djb2(char *str);
}


struct Emaths::v3
{
    union
    {
        float e[3];
        struct { float x, y, z; };
    };

    v3()                             : e{0 , 0 , 0 } {}
    v3(float e0, float e1, float e2) : e{e0, e1, e2} {}
    v3(float e0)                     : e{e0, e0, e0} {}

    inline Emaths::v3 operator-()   const              { return Emaths::v3(-e[0], -e[1], -e[2]); }
    inline float operator[] (int i) const              { return e[i]; }
    inline float &operator[](int i)                    { return e[i]; }
    inline Emaths::v3 &operator+=(const Emaths::v3 &v) { e[0] += v.e[0]; e[1] += v.e[1]; e[2] += v.e[2]; return *this; }
    inline Emaths::v3 &operator-=(const Emaths::v3 &v) { e[0] -= v.e[0]; e[1] -= v.e[1]; e[2] -= v.e[2]; return *this;}
    inline Emaths::v3 &operator*=(const float t)       { e[0] *= t; e[1] *= t; e[2] *= t; return *this; }
    inline Emaths::v3 &operator/=(const float t)       { return *this *= 1 / t; }
    inline float length_squared() const                { return e[0] * e[0] + e[1] * e[1] + e[2] * e[2]; }
    inline float length()         const                { return sqrt(length_squared()); }
    inline Emaths::v3 floor()                          { return Emaths::v3(int(x), int(y), int(z)); }
    inline Emaths::v3 normalize()                      { float L = length(); 
                                                  Emaths::v3 Res = (L > 0) ? Emaths::v3(x / L, y / L, z / L) : Emaths::v3(0, 0, 0); 
                                                  return Res; }
};

struct Emaths::iv3
{
    int x, y, z;
    iv3()                    : x(0),   y(0),   z(0){};   
    iv3(int x, int y, int z) : x(x),   y(y),   z(z){};
    inline Emaths::iv3 &operator+=(const Emaths::iv3 &v) { x += v.x; y += v.y; z += v.z; return *this; }
    inline Emaths::iv3 &operator-=(const Emaths::iv3 &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    inline Emaths::iv3 &operator*=(const int t)          { x *= t;   y *= t;   z *= t;   return *this; }
    inline Emaths::iv3 &operator/=(const int t)          { x /= t;   y /= t;   z /= t;   return *this; }
    inline int length_squared() const                    { return x * x + y * y + z * z; }
    inline int length()         const                    { return sqrt(length_squared()); }
    inline int dot(Emaths::iv3 V)                        { return x * V.x + y * V.y + z * V.z;   }
    inline Emaths::iv3 normalize()                       { int L = length(); 
                                                    Emaths::iv3 Res = (L > 0) ? Emaths::iv3(x / L, y / L, z / L) : Emaths::iv3(0, 0, 0); 
                                                    return Res; }
};
struct Emaths::iv2
{
    int x, y;
    iv2()             : x(0), y(0){};   
    iv2(int x, int y) : x(x), y(y){};
    inline Emaths::iv2 &operator+=(const Emaths::iv2 &v) { x += v.x; y += v.y; return *this; }
    inline Emaths::iv2 &operator-=(const Emaths::iv2 &v) { x -= v.x; y -= v.y; return *this; }
    inline Emaths::iv2 &operator*=(const int t)          { x *= t;   y *= t;   return *this; }
    inline Emaths::iv2 &operator/=(const int t)          { x /= t;   y /= t;   return *this; }
};

struct Emaths::v2
{
    union
    {
        float e[2];
        struct { float x, y; };
    };

    v2() : e{0, 0} {}
    v2(float e0, float e1) : e{e0, e1} {}
    v2(float e) : e{e, e} {}

    inline Emaths::v2 operator-()  const  { return Emaths::v2(-e[0], -e[1]); }
    inline float operator [](int i) const { return e[i]; }
    inline float &operator[](int i)       { return e[i]; }

    inline Emaths::v2 &operator+=(const Emaths::v2 &v) { e[0] += v.e[0]; e[1] += v.e[1]; return *this; }
    inline Emaths::v2 &operator-=(const Emaths::v2 &v) { e[0] -= v.e[0]; e[1] -= v.e[1];  return *this;}
    inline Emaths::v2 &operator*=(const float t)       { e[0] *= t;e[1] *= t; return *this;}
    inline Emaths::v2 &operator/=(const float t)       { return *this *= 1 / t;}

    inline float length_squared()     const { return e[0] * e[0] + e[1] * e[1];    }
    inline float length()             const { return sqrt(length_squared());       }
    inline float dot(Emaths::v2 V)          { return (float)(x * V.x + y * V.y);   }
    inline Emaths::v2 perpendicular()       { return Emaths::v2(y, -x);            }
    inline float perpdot(Emaths::v2 V)      { return (float)(x * V.y - y * V.x);   }
    inline float cross(Emaths::v2 V)        { return (float)(x * V.y - y * V.x);   }
    inline Emaths::v2 hadamard(Emaths::v2 V){ Emaths::v2 result(x * V.x, y * V.y); }
    inline float angle(Emaths::v2 V)        { return atan2(x * V.y - y * V.x, x * V.x + y * V.y); }// returns signed angle in radians
    inline Emaths::v2 normalize()           { float L = length(); 
                                       Emaths::v2 Res = (L > 0) ? Emaths::v2(x / L, y / L) : Emaths::v2(0, 0); 
                                       return Res; }
};

inline Emaths::v2 operator-(Emaths::v2 a, Emaths::v2 b)           { return Emaths::v2(a.x - b.x, a.y - b.y); }
inline Emaths::v2 operator+(Emaths::v2 a, Emaths::v2 b)           { return Emaths::v2(a.x + b.x, a.y + b.y); }
inline Emaths::v2 operator-(Emaths::v2 a, float b     )           { return Emaths::v2(a.x - b  , a.y - b  ); }
inline Emaths::v2 operator+(Emaths::v2 a, float b     )           { return Emaths::v2(a.x + b  , a.y + b  ); }
inline Emaths::v2 operator*(Emaths::v2 a, float b     )           { return Emaths::v2(a.x * b  , a.y * b  ); }
inline Emaths::v2 operator*(Emaths::v2 a, Emaths::v2 b)           { return Emaths::v2(a.x * b.x, a.y * b.y); }
inline Emaths::v2 operator*(float b     , Emaths::v2 a)           { return Emaths::v2(a.x * b  , a.y * b  ); }
inline Emaths::v2 operator/(Emaths::v2 a, float b     )           { return Emaths::v2(a.x / b  , a.y / b  ); }
inline Emaths::v2 operator/(Emaths::v2 a, Emaths::v2 b)           { return Emaths::v2(a.x / b.x, a.y / b.y); }
inline static float Emaths::dot    (Emaths::v2 A , Emaths::v2 B ) { return A.x * B.x + A.y * B.y;            }
inline static float Emaths::perpdot(Emaths::v2 A , Emaths::v2 B ) { return A.x * B.y - A.y * B.x;            }
inline static bool operator==      (Emaths::v2 A , Emaths::v2 B ) { return A.x == B.x && A.y == B.y;         }
inline static bool operator!=      (Emaths::v2 A , Emaths::v2 B ) { return A.x != B.x || A.y != B.y;         }
inline static Emaths::v2 rand_vector(float length) { return Emaths::v2(Emaths::rand_range(-100, 100) *0.01f *length, Emaths::rand_range(-100, 100) *0.01f *length); }

inline Emaths::v3 operator-  (Emaths::v3 a, Emaths::v3 b)    { return Emaths::v3(a.x - b.x, a.y - b.y, a.z - b.z);  }
inline Emaths::v3 operator+  (Emaths::v3 a, Emaths::v3 b)    { return Emaths::v3(a.x + b.x, a.y + b.y, a.z + b.z);  }
inline Emaths::v3 operator-  (Emaths::v3 a, float b     )    { return Emaths::v3(a.x - b  , a.y - b  , a.z - b  );  }
inline Emaths::v3 operator+  (Emaths::v3 a, float b     )    { return Emaths::v3(a.x + b  , a.y + b  , a.z + b  );  }
inline Emaths::v3 operator*  (Emaths::v3 a, float b     )    { return Emaths::v3(a.x * b  , a.y * b  , a.z * b  );  }
inline Emaths::v3 operator*  (Emaths::v3 a, int b       )    { return Emaths::v3(a.x * b  , a.y * b  , a.z * b  );  }
inline Emaths::v3 operator*  (Emaths::v3 a, Emaths::v3 b)    { return Emaths::v3(a.x * b.x, a.y * b.y, a.z * b.z);  }
inline Emaths::v3 operator*  (float b     , Emaths::v3 a)    { return Emaths::v3(a.x * b  , a.y * b  , a.z * b  );  }
inline Emaths::v3 operator/  (Emaths::v3 a, float b     )    { return Emaths::v3(a.x / b  , a.y / b  , a.z / b  );  }
inline Emaths::v3 operator/  (Emaths::v3 a, Emaths::v3 b)    { return Emaths::v3(a.x / b.x, a.y / b.y, a.z / b.z);  }
inline static bool operator==(Emaths::v3 A , Emaths::v3 B )  { return A.x == B.x && A.y == B.y && A.z == B.z;       }
inline static bool operator!=(Emaths::v3 A , Emaths::v3 B )  { return A.x != B.x || A.y != B.y || A.z != B.z;       }
inline static float Emaths::dot(Emaths::v3 A , Emaths::v3 B) { return A.x *  B.x +  A.y *  B.y +  A.z *  B.z;       }

inline Emaths::iv3 operator- (Emaths::iv3 a, Emaths::iv3 b) { return Emaths::iv3(a.x - b.x, a.y - b.y, a.z - b.z); }
inline Emaths::iv3 operator+ (Emaths::iv3 a, Emaths::iv3 b) { return Emaths::iv3(a.x + b.x, a.y + b.y, a.z + b.z); }
inline Emaths::iv3 operator- (Emaths::iv3 a, int b        ) { return Emaths::iv3(a.x - b  , a.y - b  , a.z - b  ); }
inline Emaths::iv3 operator+ (Emaths::iv3 a, int b        ) { return Emaths::iv3(a.x + b  , a.y + b  , a.z + b  ); }
inline Emaths::iv3 operator* (Emaths::iv3 a, int b        ) { return Emaths::iv3(a.x * b  , a.y * b  , a.z * b  ); }
inline Emaths::iv3 operator* (Emaths::iv3 a, Emaths::iv3 b) { return Emaths::iv3(a.x * b.x, a.y * b.y, a.z * b.z); }
inline Emaths::iv3 operator* (int b        , Emaths::iv3 a) { return Emaths::iv3(a.x * b  , a.y * b  , a.z * b  ); }
inline Emaths::iv3 operator/ (Emaths::iv3 a, int b        ) { return Emaths::iv3(a.x / b  , a.y / b  , a.z / b  ); }
inline Emaths::iv3 operator/ (Emaths::iv3 a, Emaths::iv3 b) { return Emaths::iv3(a.x / b.x, a.y / b.y, a.z / b.z); }
inline static bool operator==(Emaths::iv3 A, Emaths::iv3 B) { return A.x == B.x && A.y == B.y && A.z == B.z;}
inline static bool operator!=(Emaths::iv3 A, Emaths::iv3 B) { return A.x != B.x || A.y != B.y || A.z != B.z;}

inline Emaths::iv2 operator- (Emaths::iv2 a, Emaths::iv2 b) { return Emaths::iv2(a.x - b.x, a.y - b.y); }
inline Emaths::iv2 operator+ (Emaths::iv2 a, Emaths::iv2 b) { return Emaths::iv2(a.x + b.x, a.y + b.y); }
inline Emaths::iv2 operator- (Emaths::iv2 a, int b        ) { return Emaths::iv2(a.x - b  , a.y - b  ); }
inline Emaths::iv2 operator+ (Emaths::iv2 a, int b        ) { return Emaths::iv2(a.x + b  , a.y + b  ); }
inline Emaths::iv2 operator* (Emaths::iv2 a, int b        ) { return Emaths::iv2(a.x * b  , a.y * b  ); }
inline Emaths::iv2 operator* (Emaths::iv2 a, Emaths::iv2 b) { return Emaths::iv2(a.x * b.x, a.y * b.y); }
inline Emaths::iv2 operator* (int b        , Emaths::iv2 a) { return Emaths::iv2(a.x * b  , a.y * b  ); }
inline Emaths::iv2 operator/ (Emaths::iv2 a, int b        ) { return Emaths::iv2(a.x / b  , a.y / b  ); }
inline Emaths::iv2 operator/ (Emaths::iv2 a, Emaths::iv2 b) { return Emaths::iv2(a.x / b.x, a.y / b.y); }
inline static bool operator==(Emaths::iv2 A, Emaths::iv2 B) { return A.x == B.x && A.y == B.y  ; }
inline static bool operator!=(Emaths::iv2 A, Emaths::iv2 B) { return A.x != B.x || A.y != B.y  ; }

inline Emaths::iv2 Emaths::_iv2(Emaths::v2 A ) { Emaths::iv2 B; B.x = A.x; B.y = A.y;            return B;} 
inline Emaths::v2  Emaths::_v2 (Emaths::iv2 A) { Emaths::v2  B; B.x = A.x; B.y = A.y;            return B;}
inline Emaths::iv3 Emaths::_iv3(Emaths::v3 A ) { Emaths::iv3 B; B.x = A.x; B.y = A.y; B.z = A.z; return B;}
inline Emaths::v3  Emaths::_v3 (Emaths::iv3 A) { Emaths::v3  B; B.x = A.x; B.y = A.y; B.z = A.z; return B;}

float Emaths::clamp(float input, float min, float max)
{
    if (input < min)
        return min;
    else if (input > max)
        return max;
    else
        return input;
}
void Emaths::clamp(float *input, float min, float max)
{
    if (*input < min)
        *input = min;
    else if (*input > max)
        *input = max;
}
float Emaths::clamp(int input, int min, int max = INT_MAX)
{
    if (input < min)
        return min;
    else if (input > max)
        return max;
    else
        return input;
}
void Emaths::clamp(int *input, int min, int max = INT_MAX)
{
    if (*input < min)
        *input = min;
    else if (*input > max)
        *input = max;
}

struct Emaths::RandState
{
    uint64_t seed;
    bool initialized = false;
};
Emaths::RandState RANDSTATE;

uint64_t Emaths::rand_XOR()
{
    if (!RANDSTATE.initialized)
    {
        RANDSTATE.seed = time(NULL);
        RANDSTATE.initialized = true;
    }
    uint64_t x = RANDSTATE.seed;
    x ^= x << 9;
    x ^= x >> 5;
    x ^= x << 15;
    return RANDSTATE.seed = x;
}

int *Emaths::rand_array(int count)
{
    int* ptr = (int*)malloc(count * sizeof(int));
    srand(time(NULL));
    for (int i = 0; i < count; i++) ptr[i] = rand();
    return ptr;
}
int           Emaths::rand_range(int min, int max)          { return min + Emaths::rand_XOR() % (max - min + 1); } // Returns a random real in [min,max].
inline double Emaths::random_double()                       { return rand() / (RAND_MAX + 1.0);                  } // Returns a random real in [0,1].
inline double Emaths::random_double(double min, double max) { return min + (max - min) * random_double();        } // Returns a random real in [min,max].

float Emaths::v2_distance_2Points(Emaths::v2 A, Emaths::v2 B)
{
    return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

Emaths::v2 Emaths::unitvec_AtoB(Emaths::v2 A, Emaths::v2 B)
{
    float n = Emaths::v2_distance_2Points(A, B);
    return ((B - A) / n);
}

float Emaths::signed_angle_v2(Emaths::v2 A, Emaths::v2 B)
{
    return atan2(A.x * B.y - A.y * B.x, A.x * B.x + A.y * B.y);
}

Emaths::v2 Emaths::Rotate2D(Emaths::v2 P, float sine, float cosine)
{
    return Emaths::v2(Emaths::v2(cosine, -sine).dot(P), Emaths::v2(sine, cosine).dot(P));
}

Emaths::v2 Emaths::Rotate2D(Emaths::v2 P, float Angle)
{
    float sine = sin(Angle);
    float cosine = cos(Angle);
    return Emaths::Rotate2D(P, sine, cosine);
}
Emaths::v2 Emaths::Rotate2D(Emaths::v2 p, Emaths::v2 o, float angle)
{
    // Demonstration: https://www.desmos.com/calculator/8aaegifsba
    float s = sin(angle);
    float c = cos(angle);

    float x = (p.x - o.x) * c - (p.y - o.y) * s + o.x;
    float y = (p.x - o.x) * s + (p.y - o.y) * c + o.y;

    return Emaths::v2(x, y);
}

Emaths::v2 Emaths::Reflection2D(Emaths::v2 P, Emaths::v2 N)
{
    return P - 2 * N.dot(P) * N;
}

bool Emaths::PointInRectangle(Emaths::v2 P, Emaths::v2 A, Emaths::v2 B, Emaths::v2 C)
{
    Emaths::v2 M = P;
    Emaths::v2 AB = B - A;
    Emaths::v2 BC = C - B;
    Emaths::v2 AM = M - A;
    Emaths::v2 BM = M - B;
    if (0 <= AB.dot(AM) && AB.dot(AM) <= AB.dot(AB) &&
        0 <= BC.dot(BM) && BC.dot(BM) <= BC.dot(BC))
        return true;
    else
        return false;
}

int Emaths::sign(float x)
{   
    if (x > 0)
        return 1;
    else if (x < 0)
        return -1;
    else
        return 0;
}

static float abso(float F) { return F > 0 ? F : -F; };

uint32_t Emaths::hash_djb2(char *str)
{
    uint32_t hash = 5381;
    int c;
    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    return hash;
}

