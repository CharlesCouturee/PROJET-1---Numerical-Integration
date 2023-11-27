#include <string>
#include "Integrator.hpp"

class Midpoint : public Integrator {
public:
    std::string getName() {
        return "midpoint";
    }

    float* tmp;
    int tmplength;

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 4, implement midpoint method

        // You will probably want a temporary array in this method 
        VectorXf dpdt(n);

        // Create the temp vector
        VectorXf temp(n);

        // Setup dpdt
        derivs->derivs(t, p, dpdt);

        // Follow formula and update ststae of the system
        for (int i = 0; i < n; i++)
        {
            temp[i] = p[i] + ( (h * dpdt[i]) / 2.f);
        }

        derivs->derivs(t + h / 2.f, temp, temp);

        for (int j = 0; j < n; j++)
        {
            pout[j] = p[j] + (h * temp[j]);
        }

    }

};
