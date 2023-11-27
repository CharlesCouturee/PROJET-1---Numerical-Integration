#include <string>
#include "Integrator.hpp"

class ModifiedMidpoint : public Integrator {
public:

    std::string getName() {
        return "modified midpoint";
    }

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 5, implmement the modified midpoint (2/3) method.
        // see also efficient memory management suggestion in provided code for the Midpoint method.
        VectorXf dpdt(n);

        // Make a temp vector
        VectorXf temp(n);

        // Setup dpdt
        derivs->derivs(t, p, dpdt);

        // Apply the formula and update the state of the system
        for (int i = 0; i < n; i++)
        {
            temp[i] = p[i] + ((h * dpdt[i]) * (2.f / 3.f));
        }

        derivs->derivs(t * h * (2.f / 3.f), temp, temp);

        for (int j = 0; j < n; j++)
        {
            pout[j] = p[j] + (h * temp[j]);
        }

    }

};