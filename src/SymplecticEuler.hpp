#include <string>
#include "Integrator.hpp"


class SymplecticEuler : public Integrator {
public:
    std::string getName() {
        return "symplectic Euler";
    }

    void step(VectorXf &p, int n, float t, float h, VectorXf &pout, Function* derivs) {
        // TODO: Objective 7, complete the symplectic Euler integration method.
        // note you'll need to know how p is packed to properly implement this, so go
        // look at ParticleSystem.getPhaseSpace()
   
        VectorXf dpdt(n);

        derivs->derivs(t, p, dpdt);

        for (int i = 0; i < n; i += 4)
        {
            dpdt[i + 2] *= h;
            dpdt[i + 3] *= h;
        }

        for (int i = 0; i < n; i += 4)
        {
            pout[i + 2] = p[i + 2] + dpdt[i + 2];
            pout[i] = p[i] + pout[i + 2] * h;

            pout[i + 3] = p[i + 3] + dpdt[i + 3];
            pout[i + 1] = p[i + 1] + pout[i + 3] * h;
        }
    }

};