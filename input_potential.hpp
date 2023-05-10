#ifndef INPUT_POTENTIAL_HPP
#define INPUT_POTENTIAL_HPP
#include <cmath>
double input_potential(double x)
{
        double U_0=3;
        double L = 1;
        return U_0* pow( (pow((x/L), 2) - 1), 2);
}
#endif
