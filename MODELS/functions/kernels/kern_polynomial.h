#ifndef KERN_POLYNOMIAL_H
#define KERN_POLYNOMIAL_H
#include "kern_base.h"
/** \brief Polynomial kernel.
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \big( \gamma*\mathbf{x}^T \mathbf{y} +
 * \mathrm{c} \big)^{\mathrm{d}}
 * \f]
 *
 *  Required Config key: \ref MPARAMS <int> d  <double> \f$ \gamma \f$  c
 *
 *  @see Kern_Base
 */
class Kern_Polynomial : public virtual Kern_Base {
    public:
        Kern_Polynomial ();
        Kern_Polynomial (const Config &c);
        int d; // power
        double gamma;   // scalling
        double c;   // shift
        /**
         * Label used for this class
         */
        std::string label = "Kern_Polynomial";
        std::string get_label();
        double operator() (const aed_type2& b, const aed_type2& af);
        aed_type2 derivative(const aed_type2& b, const aed_type2& af);
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff);

};
#endif
