#include <string>
#include "Integrator.hpp"

class RK4 : public Integrator {
public:
    std::string getName() {
        return "RK4";
    }

    void step(VectorXf& p, int n, float t, float h, VectorXf& pout, Function* derivs) {
        // TODO: Objective 6, implement the RK4 integration method
        // see also efficient memory management suggestion in provided code for the Midpoint method.

        // Make necessary vectors for the formula
        VectorXf k_1(n);
        VectorXf k_2(n);
        VectorXf k_3(n);
        VectorXf k_4(n);

        // Apply formula and update system
        derivs->derivs(t, p, k_1);
        for (int i = 0; i < n; i++)
        {
            k_2[i] = p[i] + (h * k_1[i] / 2.f);
        }

        derivs->derivs(t + h / 2.f, k_2, k_2);
        for (int j = 0; j < n; j++)
        {
            k_3[j] = p[j] + (h * k_2[j] / 2.f);
        }

        derivs->derivs(t + h / 2.f, k_3, k_3);
        for (int k = 0; k < n; k++)
        {
            k_4[k] = p[k] + (h * k_3[k]);
        }

        derivs->derivs(t + h, k_4, k_4);
        for (int l = 0; l < n; l++)
        {
            pout[l] = p[l] + (h / 6.f) * (k_1[l] + (2 * k_2[l]) + (2 * k_3[l]) + k_4[l]);
        }
    }
};