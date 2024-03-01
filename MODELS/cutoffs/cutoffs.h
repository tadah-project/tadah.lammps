#ifndef CUTOFFS_H
#define CUTOFFS_H

#include <string>
#include "../../CORE/config/config.h"
#include "../../CORE/registry.h"

/** Base class to be inherited by all cutoffs */
class Cut_Base {
    public:
        virtual ~Cut_Base();
        virtual std::string label()=0;
        virtual double calc(double r)=0;
        virtual double get_rcut()=0;
        virtual void set_rcut(const double r)=0;
        virtual double get_rcut_sq()=0;
        virtual double calc_prime(double r)=0;
        void test_rcut(const double r);
};

/** \brief Dummy cutoff function
 * \f[
 * f_c(r) =
 * \begin{equation}
 * \begin{cases}
 *  1 & \text{if  } r \leq r_c\\
 *  0 & \text{otherwise}\\
 * \end{cases}
 * \end{equation}
 * \f]
 */
class Cut_Dummy : public Cut_Base {
    private:
        std::string lab = "Cut_Dummy";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Dummy();
        Cut_Dummy(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/** \brief Cos cutoff function
 *
 * \f[
 * f_c(r) =
 * \begin{equation}
 * \begin{cases}
 *  \frac{1}{2}\big[ \cos\big( \frac{\pi r}{r_c} \big)+1 \big] & \text{if  } r \leq r_c\\
 *  0 & \text{otherwise}\\
 * \end{cases}
 * \end{equation}
 * \f]
 *
 * <div class="csl-entry">Behler, J., Parrinello, M. (2007).
 * Generalized neural-network representation of high-dimensional
 * potential-energy surfaces. <i>Physical Review Letters</i>,
 * <i>98</i>(14), 146401. https://doi.org/10.1103/PhysRevLett.98.146401</div>
 */
class Cut_Cos : public Cut_Base {
    private:
        std::string lab = "Cut_Cos";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Cos();
        Cut_Cos(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/** \brief Tanh cutoff function.
 * \f[
 * f_c(r) =
 * \begin{equation}
 *  \begin{cases}
 *   \tanh^3\big( 1 -\frac{r}{r_c} \big)  & \text{if  } r \leq r_c\\
 *   0 & \text{otherwise}\\
 *  \end{cases}
 * \end{equation}
 * \f]
 *
 * <div class="csl-entry">Behler, J. (2011). Atom-centered symmetry
 * functions for constructing high-dimensional neural network potentials.
 * <i>J. Chem. Phys.</i>, <i>134</i>(7), 074106.
 * https://doi.org/10.1063/1.3553717</div>
 */
class Cut_Tanh : public Cut_Base {
    private:
        std::string lab = "Cut_Tanh";
        double rcut, rcut_sq, rcut_inv;
    public:
        Cut_Tanh();
        Cut_Tanh(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};

/** \brief Polynomial-2 cutoff function.
 *
 * \f[
 * f_c(r) =
 * \begin{equation}
 *  \begin{cases}
 *   1 & \text{if  } r \leq (r_c-1)\\
 *   r^3(r(15-6r)-10)+1  & \text{if  } (r_c-1) < r \leq r_c\\
 *   0 & \text{otherwise}\\
 *  \end{cases}
 * \end{equation}
 * \f]

 *<div class="csl-entry">Singraber, A., Rg Behler, J., Dellago, C. (2019).
 Library-Based LAMMPS Implementation of High-Dimensional
 Neural Network Potentials. https://doi.org/10.1021/acs.jctc.8b00770</div>
 */
class Cut_Poly2 : public Cut_Base {
    private:
        std::string lab = "Cut_Poly2";
        double rcut, rcut_sq, rcut_inv;
        double rcut_inner;
    public:
        Cut_Poly2();
        Cut_Poly2(double rcut);
        std::string label() ;
        void set_rcut(const double r);
        double get_rcut();
        double get_rcut_sq();
        double calc(double r);
        double calc_prime(double r);
};
//template<> inline Registry<Cut_Base,double>::Map Registry<Cut_Base,double>::registry{};
#endif
