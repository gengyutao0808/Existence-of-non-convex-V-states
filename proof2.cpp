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
const double delta2 = 0.04;
const int m2 = 20000;
const int n2 = 200000;

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
    return intervalcos(z)*(Rx*dddRy+dRx*ddRy-dddRx*Ry-ddRx*dRy) + intervalsin(z)*(Interval(2,2)*dRx*dRy+Interval(2,2)*ddRx*ddRy+ddRx*Ry+Rx*ddRy+dddRx*dRy-dRx*dddRy);
}

Interval intervaldD(const Interval& x, const Interval& z) {
    Interval Rx = intervalR(x);
    Interval Ry = intervalR(x - z);
    Interval dRx = intervaldR(x);
    Interval dRy = intervaldR(x - z);
    Interval ddRx = intervalddR(x);
    Interval ddRy = intervalddR(x - z);
    return Interval(2,2)*(Rx*ddRx+dRx*dRx+dRy*dRy+Ry*ddRy-intervalcos(z)*(ddRx*Ry+dRx*dRy)-intervalsin(z)*(dRx*dRy+Rx*ddRy));
}


int main() {
    // 开始计时
    auto start = chrono::high_resolution_clock::now();

    // 生成X2数组
    vector<double> X2;
    X2.reserve(m2 + 1);
    const double x_step = M_PI / (6 * m2);
    for (int i = 0; i <= m2; ++i) {
        X2.push_back(i * x_step);
    }

    // 生成Z2数组
    vector<double> Z2;
    Z2.reserve(n2 + 1);
    const double z_start = delta2;
    const double z_end = 2 * M_PI - delta2;
    const double z_step = (z_end - z_start) / n2;
    for (int i = 0; i <= n2; ++i) {
        Z2.push_back(z_start + i * z_step);
    }

    // 计算int2
    const double term1 = pow(delta2, 2) * (M*M4 + 2*(M*M2+M1*M1) + 3*M2*M2 + 4*M1*M3);
    const double log_term = 1 - 2 * log(2*m*delta2/M_PI);
    const double term2 = pow(delta2, 2) * (M*M3 + 2*M*M1 + 3*M1*M2);
    const double fraction1 = (M1*M2*pow(M_PI, 2) + 4*M*M1 
                           + 4*delta2*(M*M2 + pow(M1, 2))/3) / pow(m, 2);
    const double term3 = pow(delta2, 2) * (M*M2 + 2*pow(M1, 2) + M*M);
    const double fraction2 = (0.5*pow(M_PI, 2)*(M2*M2+M1*M3)+2*(M1*M1+M*M2)
                            + 2*delta2*(3*M1*M2 + M*M3)/3   ) / pow(m, 2);
    const double fraction3 = (pow((0.5*pow(M_PI, 2)*M1*M2+2*M*M1),2)+0.5*pow(delta2, 2)*pow((M*M2+M1*M1),2)
                            + 2*delta2*(0.5*pow(M_PI, 2)*M1*M2+2*M*M1)*(M*M2+M1*M1)/3  ) / pow(m, 4);
    const double int2 = (term1 * log_term + term2 * fraction1 + term3 * (fraction2+fraction3))/(4 * M_PI * o);
    
    // 输出初始信息
    cout << "delta2 = " << scientific << delta2 << endl;
    cout << "int2 = " << scientific << int2 << endl;   
    
    // 主计算参数
    const double step_size = (2*M_PI - 2*delta2)/(n2*(4 * M_PI * o));
    const Interval step(step_size, step_size);


    // 主循环
    double ddE0 = 1;
    try {
        omp_set_num_threads(56);
        #pragma omp parallel for 
        for (int l = 0; l < m2; ++l){
            // 初始化区间
            Interval x(X2[l], X2[l+1]);
            // 计算初始误差
            Interval dR_x = intervaldR(x);
            Interval ddR_x = intervalddR(x);
            Interval dddR_x = intervaldddR(x);
            Interval R_x = intervalR(x);
            Interval err = Interval(3,3)*dR_x * ddR_x + R_x * dddR_x;

            for (int k = 0; k < n2; ++k) {
                    Interval z(Z2[k], Z2[k+1]);
                    // 计算各项
                    Interval A = intervalA(x, z);
                    Interval B = intervalB(x, z);
                    Interval C = intervalC(x, z);
                    Interval D = intervalD(x, z);
                    Interval dC = intervaldC(x, z);
                    Interval dD = intervaldD(x, z);
                    // 计算当前项
                    Interval ln_A = intervalln(A);
                    Interval term1 = ln_A * dC;
                    Interval term2 = (D * C) / A;
                    Interval term3 = (B * dD) / A;
                    Interval term4 = (B * D * D) / (A*A);
                    Interval delta = term1 + term2 + term2 + term3 - term4;       
                    err = err - delta * step;     
            }
            
            // 最终输出
            // cout << l << ":" << err << endl;
            if (err.upper >= 50-int2){
                cout << x << ":upper error!" << err << endl;
                ddE0=0;
            }
            if (err.lower <= -50+int2){
                cout << x << ":lower error!" << err << endl;
                ddE0=0;
            }
        }
    } catch (const exception& e) {
        cerr << "\nError occurred: " << e.what() << endl;
        return 1;
    }

    if(ddE0 == 1){
        cout << "|ddE[x]| is bounded in 50." << endl;
    }
    if(ddE0 == 0){
        cout << "Proof is error!" << endl;
    }

    //结束计时
    auto end = chrono::high_resolution_clock::now();
    // 计算时间差并转换为秒
    chrono::duration<double> duration = end - start;
    cout << "程序运行时间: " << duration.count() << " 秒" << endl;
return 0;
}


