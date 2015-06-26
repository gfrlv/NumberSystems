#include <iostream>
#include "support.h"
#include "complex.h"
#include "quaternion.h"
#include <cmath>

typedef Quaternion<double> Quat;

Quat rotate(Quat r, double angle, Quat axis){
    Quat w = cos(0.5 * angle) * Quat(1, 0, 0, 0)
            + sin(0.5 * angle) * axis / axis.norm();
    return w * r * w.inverse();
}

int main(){
    Quat r(0, 0, 1, 0);
    double theta = -M_PI * 0.25;
    std::cout << rotate(r, theta, Quat(0, 1, 0, 0)) << std::endl;

}
