#include <bits/stdc++.h>
#include "omp.h"
#include <chrono>

using namespace std;

const int N = 30;
const double o = 1537.0 / 3750.0;
const vector<double> c = {
    1.0, 0.06990356205479509, 0.0148148550964033, 0.00456258086735692,
    0.0016620666483638713, 0.0006660059943472454, 0.0002837411049042706,
    0.00012614991045596885, 5.7871895865025096e-05, 2.7193965356809378e-05,
    1.3022964203014303e-05, 6.333013816999933e-06, 3.118981248161905e-06,
    1.552511163626916e-06, 7.798097208979299e-07, 3.947485654774547e-07,
    2.011788504215774e-07, 1.0313790964680969e-07, 5.3150561947961317e-08,
    2.7513294912126838e-08, 1.429886987865343e-08, 7.458296874209281e-09,
    3.9018181293119985e-09, 2.046434382043023e-09, 1.0758700103758797e-09,
    5.667181714328583e-10, 2.9890572119286554e-10, 1.5796673367306458e-10,
    8.330084461317743e-11, 4.4460251640104796e-11, 2.1976876953775598e-11
};

const double m = 0.941;
const double M = 1.0925;
const double M1 = 0.52;
const double M2 = 8.7;
const double M3 = 106.0;
const double M4 = 4000.0;
const double delta0 = 1.4e-5;
const double delta1 = 2e-3;
const int m0 = 100000;
const int n0 = 600000;
const int n1 = 200000;

class Interval {
public:
    double lower;
    double upper;

    Interval(double l, double u) : lower(l), upper(u) {}
    
    Interval operator+(const Interval& other) const {
        return Interval(lower + other.lower, upper + other.upper);
    }
    
    Interval operator-(const Interval& other) const {
        return Interval(lower - other.upper, upper - other.lower);
    }
    
    Interval operator*(const Interval& other) const {
        vector<double> products = {
            lower * other.lower,
            lower * other.upper,
            upper * other.lower,
            upper * other.upper
        };

        auto minmax = minmax_element(products.begin(), products.end());
        return Interval(*minmax.first, *minmax.second);
    }
    
    Interval operator/(const Interval& other) const {
        if (other.lower <= 0 && other.upper >= 0) {
            throw runtime_error("Division by interval containing zero");
        }
        vector<double> reciprocals = {1.0 / other.upper, 1.0 / other.lower};
        auto minmax = minmax_element(reciprocals.begin(), reciprocals.end());
        return *this * Interval(*minmax.first, *minmax.second);
    }

    friend ostream& operator<<(ostream& os, const Interval& ival) {
        os << "[" << fixed << setprecision(16) <<  ival.lower << ", " << fixed << setprecision(16) <<  ival.upper << "]";
        return os;
    }
};

Interval intervalR(const Interval& x) {
    Interval s(0, 0);
    for(int k=0; k<=N; ++k) {
        double theta_min = 6 * k * x.lower;
        double theta_max = 6 * k * x.upper;
        double delta_theta = theta_max - theta_min;

        if(delta_theta >= 2*M_PI) {
            s = s + Interval(-c[k], c[k]);
        } else {
            // 检查极值点
            int n_min = ceil(theta_min / (2*M_PI));
            int n_max = floor(theta_max / (2*M_PI));
            bool has_max = (n_min <= n_max);

            int m_min = ceil((theta_min - M_PI) / (2*M_PI));
            int m_max = floor((theta_max - M_PI) / (2*M_PI)); 
            bool has_min = (m_min <= m_max);

            // 端点值
            double cos_min = cos(theta_min);
            double cos_max = cos(theta_max);
            double end_min = min(cos_min, cos_max);
            double end_max = max(cos_min, cos_max);

            // 确定值域
            double min_val = has_min ? -1.0 : end_min;
            double max_val = has_max ? 1.0 : end_max;
            s = s + Interval(c[k]*min_val, c[k]*max_val);
        }
    }
    return s;
}

Interval intervaldR(const Interval& x) {
    Interval s(0, 0);
    for(int k=0; k<=N; ++k) {
        double theta_min = 6 * k * x.lower;
        double theta_max = 6 * k * x.upper;
        double delta_theta = theta_max - theta_min;

        if(delta_theta >= 2*M_PI) {
            s = s - Interval(-6*k*c[k], 6*k*c[k]);
        } else {
            // 检查极值点
            int p_min = ceil((theta_min - M_PI_2) / (2*M_PI));
            int p_max = floor((theta_max - M_PI_2) / (2*M_PI));
            bool has_max = (p_min <= p_max);

            int q_min = ceil((theta_min - 3*M_PI_2) / (2*M_PI));
            int q_max = floor((theta_max - 3*M_PI_2) / (2*M_PI));
            bool has_min = (q_min <= q_max);

            // 端点值
            double sin_min = sin(theta_min);
            double sin_max = sin(theta_max);
            double end_min = min(sin_min, sin_max);
            double end_max = max(sin_min, sin_max);

            // 确定值域
            double min_val = has_min ? -1.0 : end_min;
            double max_val = has_max ? 1.0 : end_max;
            s = s - Interval(6*k*c[k]*min_val, 6*k*c[k]*max_val);
        }
    }
    return s;
}

Interval intervalddR(const Interval& x) {
    Interval s(0, 0);
    for(int k=0; k<=N; ++k) {
        double theta_min = 6 * k * x.lower;
        double theta_max = 6 * k * x.upper;
        double delta_theta = theta_max - theta_min;

        if(delta_theta >= 2*M_PI) {
            s = s - Interval(-36*k*k*c[k], 36*k*k*c[k]);
        } else {
            // 检查极值点
            int n_min = ceil(theta_min / (2*M_PI));
            int n_max = floor(theta_max / (2*M_PI));
            bool has_max = (n_min <= n_max);

            int m_min = ceil((theta_min - M_PI) / (2*M_PI));
            int m_max = floor((theta_max - M_PI) / (2*M_PI));
            bool has_min = (m_min <= m_max);

            // 端点值
            double cos_min = cos(theta_min);
            double cos_max = cos(theta_max);
            double end_min = min(cos_min, cos_max);
            double end_max = max(cos_min, cos_max);

            // 确定值域
            double min_val = has_min ? -1.0 : end_min;
            double max_val = has_max ? 1.0 : end_max;
            s = s - Interval(36*k*k*c[k]*min_val, 36*k*k*c[k]*max_val);
        }
    }
    return s;
}

Interval intervaldddR(const Interval& x) {
    Interval s(0, 0);
    for(int k=0; k<=N; ++k) {
        double theta_min = 6 * k * x.lower;
        double theta_max = 6 * k * x.upper;
        double delta_theta = theta_max - theta_min;

        if(delta_theta >= 2*M_PI) {
            s = s + Interval(-216*k*k*k*c[k], 216*k*k*k*c[k]);
        } else {
            // 检查极值点
            int p_min = ceil((theta_min - M_PI_2) / (2*M_PI));
            int p_max = floor((theta_max - M_PI_2) / (2*M_PI));
            bool has_max = (p_min <= p_max);

            int q_min = ceil((theta_min - 3*M_PI_2) / (2*M_PI));
            int q_max = floor((theta_max - 3*M_PI_2) / (2*M_PI));
            bool has_min = (q_min <= q_max);

            // 端点值
            double sin_min = sin(theta_min);
            double sin_max = sin(theta_max);
            double end_min = min(sin_min, sin_max);
            double end_max = max(sin_min, sin_max);

            // 确定值域
            double min_val = has_min ? -1.0 : end_min;
            double max_val = has_max ? 1.0 : end_max;
            s = s + Interval(216*k*k*k*c[k]*min_val, 216*k*k*k*c[k]*max_val);
        }
    }
    return s;
}

Interval intervalsin(const Interval& x) {
    double delta = x.upper - x.lower;
    if(delta >= 2*M_PI) {
        return Interval(-1, 1);
    }

    // 检查极值点
    int p_min = ceil((x.lower - M_PI_2) / (2*M_PI));
    int p_max = floor((x.upper - M_PI_2) / (2*M_PI));
    bool has_max = (p_min <= p_max);

    int q_min = ceil((x.lower - 3*M_PI_2) / (2*M_PI));
    int q_max = floor((x.upper - 3*M_PI_2) / (2*M_PI));
    bool has_min = (q_min <= q_max);

    // 端点值
    double sin_low = sin(x.lower);
    double sin_high = sin(x.upper);
    double end_min = min(sin_low, sin_high);
    double end_max = max(sin_low, sin_high);

    return Interval(has_min ? -1.0 : end_min, 
                   has_max ? 1.0 : end_max);
}

Interval intervalcos(const Interval& x) {
    double delta = x.upper - x.lower;
    if(delta >= 2*M_PI) {
        return Interval(-1, 1);
    }

    // 检查极值点
    int n_min = ceil(x.lower / (2*M_PI));
    int n_max = floor(x.upper / (2*M_PI));
    bool has_max = (n_min <= n_max);

    int m_min = ceil((x.lower - M_PI) / (2*M_PI));
    int m_max = floor((x.upper - M_PI) / (2*M_PI));
    bool has_min = (m_min <= m_max);

    // 端点值
    double cos_low = cos(x.lower);
    double cos_high = cos(x.upper);
    double end_min = min(cos_low, cos_high);
    double end_max = max(cos_low, cos_high);

    return Interval(has_min ? -1.0 : end_min,
                   has_max ? 1.0 : end_max);
}

Interval intervalln(const Interval& x) {
    if(x.lower <= 0) {
        throw runtime_error("Logarithm of non-positive interval");
    }
    return Interval(log(x.lower), log(x.upper));
}

// 辅助组合函数
Interval intervalA(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    return Rx*Rx + Ry*Ry - Interval(2,2)*Rx*Ry*intervalcos(z);
}

Interval intervalB(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    return intervalcos(z)*(Rx*dRy - dRx*Ry) + intervalsin(z)*(Rx*Ry + dRx*dRy);
}

Interval intervalC(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    return intervalcos(z)*(Rx*ddRy - ddRx*Ry) + intervalsin(z)*(dRx*Ry + Rx*dRy + ddRx*dRy + dRx*ddRy);
}

Interval intervalD(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    return Interval(2,2)*(Rx*dRx + Ry*dRy - dRx*Ry*intervalcos(z) - Rx*dRy*intervalcos(z));
}

Interval intervalE(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    return intervalcos(z)*(dRx*dRy + dRx*dRy - Rx*ddRy + Rx*Ry) + intervalsin(z)*( dRx*Ry - dRx*ddRy - Rx*dRy - Rx*dRy);
}

Interval intervalF(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    return Interval(2,2)*(Rx*Ry*intervalsin(z) + Rx*dRy*intervalcos(z) - dRy*Ry);
}

Interval intervaldC(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    Interval dddRx = intervaldddR(x);
    Interval dddRy = intervaldddR(x - z);
    return intervalcos(z)*(Rx*dRy-Rx*dddRy+dRx*Ry+dRx*ddRy+Interval(2,2)*ddRx*dRy) - intervalsin(z)*(Interval(2,2)*Rx*ddRy +dRx*dRy+dRx*dddRy+ddRx*ddRy-ddRx*Ry);
}

Interval intervaldD(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    return Interval(2,2)*(intervalcos(z)*(dRx*dRy+Rx*ddRy)+intervalsin(z)*(dRx*Ry+Rx*dRy)-dRy*dRy-Ry*ddRy);
}

Interval intervaldE(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    Interval dddRx = intervaldddR(x);
    Interval dddRy = intervaldddR(x - z);
    return intervalsin(z)*(Interval(3,3)*Rx*ddRy - Interval(3,3)*dRx*dRy + dRx*dddRy - Rx*Ry) 
         + intervalcos(z)*(dRx*Ry + Rx*dddRy - Interval(3,3)*dRx*ddRy - Interval(3,3)*Rx*dRy);
}

Interval intervaldF(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    return Interval(2,2)*(dRy*dRy + Ry*ddRy + Rx*(intervalcos(z)*(Ry-ddRy)-intervalsin(z)*(dRy+ddRy)));
}

double R(const double& x) {
    double s = 0.0;
    for (int i = 0; i <= N; ++i) {
        s += c[i] * cos(6 * i * x);
    }
    return s;
}

double dR(const double& x) {
    double s = 0.0;
    for (int i = 0; i <= N; ++i) {
        s += -6 * i * c[i] * sin(6 * i * x);
    }
    return s;
}

double ddR(const double& x) {
    double s = 0.0;
    for (int i = 0; i <= N; ++i) {
        s += -36 * i * i * c[i] * cos(6 * i * x);
    }
    return s;
}

double FF(const double& x, const double& z) {
    double Rx = 0.0;
    double Ry = 0.0;
    for (int i = 0; i <= N; ++i) {
        Rx += c[i] * cos(6 * i * x);
        Ry += c[i] * cos(6 * i * (x - z));
    }

    double dRx = 0.0;
    double dRy = 0.0;
    for (int i = 0; i <= N; ++i) {
        dRx += -6 * i * c[i] * sin(6 * i * x);
        dRy += -6 * i * c[i] * sin(6 * i * (x - z));
    }

    double f0 = Rx * Rx + Ry * Ry - 2 * Rx * Ry * cos(z);
    double f1 = cos(z) * (Rx * dRy - dRx * Ry) + sin(z) * (Rx * Ry + dRx * dRy);

    double F1 = log(f0) * f1 ;
    return F1;
}

double dFF(const double& x, const double& z) {
    double Rx = 0.0;
    double Ry = 0.0;
    for (int i = 0; i <= N; ++i) {
        Rx += c[i] * cos(6 * i * x);
        Ry += c[i] * cos(6 * i * (x - z));
    }

    double dRx = 0.0;
    double dRy = 0.0;
    for (int i = 0; i <= N; ++i) {
        dRx += -6 * i * c[i] * sin(6 * i * x);
        dRy += -6 * i * c[i] * sin(6 * i * (x - z));
    }

    double ddRx = 0.0;
    double ddRy = 0.0;
    for (int i = 0; i <= N; ++i) {
        ddRx += -36 * i * i * c[i] * cos(6 * i * x);
        ddRy += -36 * i * i * c[i] * cos(6 * i * (x - z));
    }

    double f0 = Rx * Rx + Ry * Ry - 2 * Rx * Ry * cos(z);
    double f1 = cos(z) * (Rx * dRy - dRx * Ry) + sin(z) * (Rx * Ry + dRx * dRy);
    double f2 = cos(z) * (Rx*ddRy-ddRx*Ry) + sin(z) * (dRx*Ry+Rx*dRy+ddRx*dRy+dRx*ddRy);
    double f3 = 2 * (Rx*dRx+Ry*dRy-dRx*Ry*cos(z)-Rx*dRy*cos(z));

    double F1 = log(f0)*f2+f1*f3/f0 ;
    return F1;
}

int main() {
    // 开始计时
    auto start = chrono::high_resolution_clock::now();

    // 生成X0数组
    vector<double> X0;
    X0.reserve(m0 + 1);
    const double x0_step = M_PI / (6 * m0);
    for (int i = 0; i <= m0; ++i) {
        X0[i] = i * x0_step ;
    }

    // 生成Z0数组
    vector<double> Z0;
    Z0.reserve(n0 + 1);
    const double z0_start = delta0 + 0.13335;
    const double z0_end = 2 * M_PI - delta0 - 0.13335;
    const double z0_step = (z0_end - z0_start) / 60000;
    for (int i = 0; i <= 40000; ++i) {
        Z0[i] = delta0 + i * 1e-10;
        Z0[n0-i] = 2 * M_PI - delta0 - i * 1e-10;
    }
    for (int i = 0; i <= 30000; ++i) {
        Z0[40000+i] = delta0 + 4e-6 + i * 2e-10;
        Z0[n0-40000-i] = 2 * M_PI - delta0 - 4e-6 - i * 2e-10;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[70000+i] = delta0 + 1e-5 + i * 3e-10;
        Z0[n0-70000-i] = 2 * M_PI - delta0 - 1e-5 - i * 3e-10;
    }
    for (int i = 0; i <= 10000; ++i) {
        Z0[90000+i] = delta0 + 0.000016 + i * 4e-10;
        Z0[n0-90000-i] = 2 * M_PI - delta0 - 0.000016 - i * 4e-10;
    }

    for (int i = 0; i <= 20000; ++i) {
        Z0[100000+i] = delta0 + 0.00002 + i * 6e-10;
        Z0[n0-100000-i] = 2 * M_PI - delta0 - 0.00002 - i * 6e-10;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[120000+i] = delta0 + 0.000032 + i * 9e-10;
        Z0[n0-120000-i] = 2 * M_PI - delta0 - 0.000032 - i * 9e-10;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[140000+i] = delta0 + 0.00005 + i * 3e-9;
        Z0[n0-140000-i] = 2 * M_PI - delta0 - 0.00005 - i * 3e-9;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[160000+i] = delta0 + 0.00011 + i * 1.2e-8;
        Z0[n0-160000-i] = 2 * M_PI - delta0 - 0.00011 - i * 1.2e-8;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[180000+i] = delta0 + 0.00035 + i * 5e-8;
        Z0[n0-180000-i] = 2 * M_PI - delta0 - 0.00035 - i * 5e-8;
    }

    for (int i = 0; i <= 20000; ++i) {
        Z0[200000+i] = delta0 + 0.00135 + i * 1e-7;
        Z0[n0-200000-i] = 2 * M_PI - delta0 - 0.00135 - i * 1e-7;
    }
    for (int i = 0; i <= 30000; ++i) {
        Z0[220000+i] = delta0 + 0.00335 + i * 1e-6;
        Z0[n0-220000-i] = 2 * M_PI - delta0 - 0.00335 - i * 1e-6;
    }
    for (int i = 0; i <= 20000; ++i) {
        Z0[250000+i] = delta0 + 0.03335 + i * 5e-6;
        Z0[n0-250000-i] = 2 * M_PI - delta0 - 0.03335 - i * 5e-6;
    }

    for (int i = 0; i <=60000; ++i) {
        Z0[270000+i] = z0_start + i * z0_step;
    }

    // 生成Z1数组
    vector<double> Z1;
    Z1.reserve(n1 + 1);
    const double z1_start = delta1 + 0.101;
    const double z1_end = 2 * M_PI - delta1 - 0.101;
    const double z1_step = (z1_end - z1_start) / (n1/2);
    for(int i = 0; i<=10000; ++i){
        Z1[i] = delta1 + i * 1e-7;
        Z1[n1-i] = 2 * M_PI - delta1 - i * 1e-7;
    }
    for (int i = 0; i <=40000; ++i){
        Z1[10000+i] = delta1 + 0.001 + i * 2.5e-6;
        Z1[n1-10000-i] = 2 * M_PI - delta1 - 0.001 - i * 2.5e-6;
    }
    for (int i = 0; i <=100000; ++i) {
        Z1[50000+i] = z1_start + i * z1_step;
    }

    // 计算int0
    const double term = pow(delta0, 2) * (M*M2 + 2*pow(M1, 2) + M*M);
    const double log_term = 1 - 2 * log(2*m*delta0/M_PI);
   
    const double int0 = term * log_term/(4 * M_PI * o);
    
    // 计算int1
    const double term1 = pow(delta1, 2) * (M*M3 + 2*M*M1 + 3*M1*M2);
    const double log_term1 = 1 - 2 * log(2*m*delta1/M_PI);
    const double term2 = pow(delta1, 2) * (M*M2 + 2*pow(M1, 2) + M*M);
    const double fraction = (0.5*M1*M2*pow(M_PI, 2) + 2*M*M1 
                           + 2*delta1*(M*M2 + pow(M1, 2))/3) / pow(m, 2);
    
    const double int1 = (term1 * log_term1 + term2 * fraction)/(4 * M_PI * o);
    
    // 输出初始信息
    cout << "delta0 = " << scientific << delta0 << endl;
    cout << "int0 = " << scientific << int0 << endl;
    cout << "delta1 = " << scientific << delta1 << endl;
    cout << "int1 = " << scientific << int1 << endl;
    cout << "int1*x0_step/2=" << scientific << int1*x0_step/2 << endl;

    // 主计算参数
    const Interval dEstep(-x0_step/2, x0_step/2);
    const Interval ddEstep(-25 * (x0_step/2) * (x0_step/2), 25 * (x0_step/2) * (x0_step/2));


    // 主循环
    double E0 = 1; 
    omp_set_num_threads(56);
    #pragma omp parallel for 
    for (int l = 0; l < m0; ++l){
        try{
            // 初始化区间
            double x=(X0[l]+X0[l+1])/2;
            Interval X(x,x);
            // 计算初始误差
            double dR_x = dR(x);
            double ddR_x = ddR(x);
            double R_x = R(x);
            double err0 = R_x * dR_x ;
            double err1 = dR_x * dR_x + R_x * ddR_x;
            Interval Err0(err0,err0);
            double upp0 = 0;
            double low0 = 0;
            Interval Err1(err1,err1);
            double upp1 = 0;
            double low1 = 0;
            Interval Err(0,0);

            for (int k = 0; k < n0; ++k) {
                Interval z(Z0[k], Z0[k+1]);
                
                // 计算各项
                Interval A = intervalA(X, z);
                Interval B = intervalB(X, z);
                Interval E = intervalE(X, z);
                Interval F = intervalF(X, z);
                Interval dE = intervaldE(X, z);
                Interval dF = intervaldF(X, z);
                
                // 计算当前项
                Interval ln_A = intervalln(A);
                Interval term1 = ln_A * dE;
                Interval term2 = Interval(2,2)*(E * F) / A;
                Interval term3 = B*(A*dF-F*F) / (A*A);
                Interval delta = term1 + term2 + term3;
                double ff = FF(x, (Z0[k]+Z0[k+1])/2);
                // 更新误差
                upp0 = upp0 + (ff*(Z0[k+1]-Z0[k]) + delta.upper*pow((Z0[k+1]-Z0[k]),3)/24)/(4 * M_PI * o);
                low0 = low0 + (ff*(Z0[k+1]-Z0[k]) + delta.lower*pow((Z0[k+1]-Z0[k]),3)/24)/(4 * M_PI * o);
            }
            Err0 = Err0 - Interval(low0,upp0) + Interval(-int0,int0);   

            for (int k = 0; k < n1; ++k) {                
                Interval z(Z1[k], Z1[k+1]);

                // 计算各项
                Interval A = intervalA(X, z);
                Interval B = intervalB(X, z);
                Interval C = intervalC(X, z);
                Interval D = intervalD(X, z);
                Interval E = intervalE(X, z);
                Interval F = intervalF(X, z);
                Interval dC = intervaldC(X, z);
                Interval dD = intervaldD(X, z);
                // 计算当前项
                Interval ln_A = intervalln(A);
                Interval term1 = ln_A * dC;
                Interval term2 = (F * C) / A;
                Interval term3 = (D*E+B*dD) / A;
                Interval term4 = (F * B * D) / (A*A);
                Interval delta1 = (term1 + term2 + term3 - term4) * Interval(-(Z1[k+1]-Z1[k])/2,(Z1[k+1]-Z1[k])/2);
                double delta2 = dFF(x,(Z1[k]+Z1[k+1])/2);

                upp1 = upp1 + (delta2 + delta1.upper)*(Z1[k+1]-Z1[k])/(4 * M_PI * o);
                low1 = low1 + (delta2 + delta1.lower)*(Z1[k+1]-Z1[k])/(4 * M_PI * o); 
            }
            Err1 = Err1 - Interval(low1,upp1) + Interval(-int1,int1);

            Err = Err0 + Err1 * dEstep + ddEstep;

            // 最终输出
            //cout << l << ":" << Err << endl;
            if (Err.upper > 3e-8){
                cout << x << ":upper error!" << Err << endl;
                E0 = 0;
            }
            if (Err.lower < -3e-8){
                cout << x << ":lower error!" << Err << endl;
                E0 = 0;
            }

        }catch (const exception& e) {
            cerr << "Thread " << omp_get_thread_num() << " error: " << e.what() << endl;            
        }        
                  
    }

    if(E0 == 1){
        cout << "|E[x]| is bounded in 3e-8." << endl;
    }
    if(E0 == 0){
        cout << "Proof is error!" << endl;
    }

    //结束计时
    auto end = chrono::high_resolution_clock::now();
    // 计算时间差并转换为秒
    chrono::duration<double> duration = end - start;
    cout << "程序运行时间: " << duration.count() << " 秒" << endl;

    return 0;
}
