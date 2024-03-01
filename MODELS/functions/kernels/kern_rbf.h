#ifndef KERN_RBF_H
#define KERN_RBF_H
#include "kern_base.h"
/** \brief Radial Basis Function kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \exp\Big( -\frac{||\mathbf{x}-\mathbf{y}||^2}
 * {2\sigma^2} \Big)
 * \f]
 *
 *  Required Config key: \ref MPARAMS <double> \f$ \sigma \f$
 *
 *  @see Kern_Base Kern_Linear Kern_Quadratic
 */
class Kern_RBF : public virtual Kern_Base {
    public:
        double sigma;
        /**
         * Label used for this class
         */
        Kern_RBF ();
        Kern_RBF (const Config &c);
        std::string label = "Kern_RBF";
        double operator() (const aed_type2& b, const aed_type2& af);
        aed_type2 derivative(const aed_type2& b, const aed_type2& af);
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff);
        std::string get_label();

    private:
        double gamma;

};
#endif
