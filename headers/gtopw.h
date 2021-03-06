#ifndef GTOPW_H
#define GTOPW_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cfenv>
#include <vector>
#include <complex>
#include <tuple>
#include <cstring>

extern "C"
{
#include <quadmath.h>
}

/* type definitions */

using usint = unsigned short int;
using uint  = unsigned int;
using llint = long long int;
using cdouble = std::complex<double>;
using quad = __float128;

constexpr int bas_lmax = 8;
constexpr double fourpi = 4.0 * M_PI;
constexpr double sqrt_4pi = sqrt(fourpi);
constexpr cdouble I(0.0, 1.0);
constexpr bool gamess_order = true;
constexpr double pi52 = pow(M_PI, 2.5);

/* loop-up tables */

constexpr int sph_siz[11] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19};
constexpr int crt_siz[11] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55};

extern int belt[500];

extern double fact[161];
extern double dfact[201];
extern double binom[101][101];
extern double omega[33][33];

extern double clmr[bas_lmax + 1][2 * bas_lmax + 1][crt_siz[bas_lmax]];
extern double xyz_norm[bas_lmax + 1][crt_siz[bas_lmax]];

/* function prototypes */

template <typename T, typename U>
inline std::complex<T> operator*(std::complex<T> lhs, const U &rhs) { return lhs *= rhs; };

inline void PrintCMatrix(const cdouble *m, const int n)
{
    /* printing function for complex packed matrix */

    std::cout << " Real component:" << std::endl;
    std::cout << std::showpos;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << real(m[i * n + j]) << " ";
        std::cout << std::endl;
    };
    std::cout << std::endl;
    std::cout << " Imag component:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << imag(m[i * n + j]) << " ";
        std::cout << std::endl;
    };
    std::cout << std::noshowpos;
    /* */
};

inline void PrintCMatrixRect(const cdouble *m, const int k, const int n)
{
    /* printing function for complex packed matrix */

    std::cout << " Real component:" << std::endl;
    std::cout << std::showpos;
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << real(m[i * k + j]) << " ";
        std::cout << std::endl;
    };
    std::cout << std::endl;
    std::cout << " Imag component:" << std::endl;
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << imag(m[i * k + j]) << " ";
        std::cout << std::endl;
    };
    std::cout << std::noshowpos;
    /* */
};

inline void PrintMatrix(const double **m, const int n)
{
    /* printing function for 2D matrix */
    std::cout << std::showpos;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << m[i][j] << " ";
        std::cout << std::endl;
    };
    std::cout << std::noshowpos;
};

inline void PrintMatrixPacked(const double *m, const int n)
{
    /* printing function for packed matrix */
    std::cout << std::showpos;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << std::setw(25) << m[i * n + j] << " ";
        std::cout << std::endl;
    };
    std::cout << std::noshowpos;
};

inline void PrintQuadMath(quad input)
{
    /* simply print the quadruple precision number */
    char buf[128];
    quadmath_snprintf(buf, sizeof buf, "%+-#*.32Qe", 32, input);
    printf("%s\n", buf);
    return;
};

void Init();

/* class prototypes */

class Keywords
{
  public:
    bool spher;
    bool cart;

    bool stvh;
    bool dip;
    bool quad;
    bool velo;
    bool eri;
    bool proj;
    bool kinc;
    bool capi;

    int write_point;
    std::vector<std::tuple<double, double, double>> points;

    std::string file1E_crt;
    std::string file1E_sph;
    std::string file2E;
    std::string proj1E;
    std::string norm1E;
    std::string path;

    double thrsh;
    double cap_on;
    double cap_mul;
    int cap_pow;

    usint n_rad;
    usint lproj_max;
    double rad_str;
    double rad_stp;

    Keywords()
    {
        stvh = false;
        dip = false;
        quad = false;
        velo = false;
        eri = false;
        proj = false;
        kinc = false;
        capi = false;

        write_point = 0;

        spher = false;
        cart = true;

        thrsh = 1.0e-14;

        n_rad = 0;
        lproj_max = -1;

        rad_str = 0.0;
        rad_stp = 0.0;
        cap_on = 10.0;
        cap_mul = 1.0;
        cap_pow = 2;
    };

    void name_me(const std::string name)
    {
        const auto file = name.substr(0, name.length() - 4);
        file1E_crt = path + "file1E_" + file + "_crt.F";
        file1E_sph = path + "file1E_" + file + "_sph.F";
        file2E = path + "file2E_" + file + ".F";
        proj1E = path + "proj1E_" + file;
        norm1E = path + "norm1E_" + file + ".F";

        std::cout << " Names of the integral files: " << std::endl;
        std::cout << "  " << file1E_crt << std::endl;
        std::cout << "  " << file1E_sph << std::endl;
        std::cout << "  " << file2E << std::endl;
        std::cout << "  " << norm1E << std::endl;
        std::cout << " (automatic assignment)" << std::endl;
    };

    ~Keywords(){};
};

template <typename T>
class arr1d
{
  public:
    int num;
    T *v;

    T &operator[](std::size_t i) { return v[i]; };
    arr1d(const int num_) : num(num_)
    {
        v = new T[num];
        std::memset(v, 0, num * sizeof(T));
    };
    ~arr1d() { delete[] v; };
    void zero() { std::memset(v, 0, num * sizeof(T)); };
};

template <typename T>
class arr2d
{
  public:
    int num_x;
    int num_y;

    T **v;

    arr2d(const int num_x_, const int num_y_) : num_x(num_x_), num_y(num_y_)
    {
        v = new T *[num_x];
        for (int i = 0; i < num_x; i++)
        {
            v[i] = new T[num_y];
            std::memset(v[i], 0, num_y * sizeof(T));
        };
    };

    ~arr2d()
    {
        for (int i = 0; i < num_x; i++)
            delete[] v[i];
        delete[] v;
    };

    void zero()
    {
        for (int i = 0; i < num_x; i++)
            std::memset(v[i], 0, num_y * sizeof(T));
    };
};

class GaussC
{
  public:
    double Ax, Ay, Az;
    double kx, ky, kz;

    int lA;
    int clen;

    std::vector<double> alphaA;
    std::vector<double> dA_re;
    std::vector<double> dA_im;
};

class Nuclei
{
  public:
    double Cx, Cy, Cz;
    double chrg;
    std::string name;
};

inline void print_math_errors(std::ostream& os) {
    if (std::fetestexcept(FE_DIVBYZERO))
        os << " FE_DIVBYZERO (pole error) reported\n";

    if (std::fetestexcept(FE_OVERFLOW))
        os << " FE_OVERFLOW (overflow error) reported\n";

    if (std::fetestexcept(FE_UNDERFLOW))
        os << " FE_UNDERFLOW (underflow error) reported\n";

    if (std::fetestexcept(FE_INEXACT))
        os << " FE_INEXACT (inexact error) reported\n";

    if (std::fetestexcept(FE_INVALID))
        os << " FE_INVALID (domain error) reported\n";
}

inline void print_if_math_errors_set(std::ostream& os) {
    os << " MATH_ERRNO is      "
        << (math_errhandling & MATH_ERRNO ? "set" : "not set") << '\n'
        << " MATH_ERREXCEPT is  "
        << (math_errhandling & MATH_ERREXCEPT ? "set" : "not set") << '\n';
}

inline void reset_math_errors() {
    std::feclearexcept(FE_ALL_EXCEPT);
}

#endif //GTOPW_H
