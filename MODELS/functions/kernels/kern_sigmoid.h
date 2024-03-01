#ifndef KERN_SIGMOID_H
#define KERN_SIGMOID_H
#include "kern_base.h"
/** \brief Sigmoid kernel.
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \tanh\big( \gamma*\mathbf{x}^T \mathbf{y}
 * + \mathrm{c} \big)
 * \f]
 *
 *  Required Config key: \ref MPARAMS <double> \f$ \gamma \f$  c
 *
 *  @see Kern_Base
 */
class Kern_Sigmoid : public virtual Kern_Base {
    public:
        double gamma;   // scalling
        double c;   // shift
        /**
         * Label used for this class
         */
        Kern_Sigmoid ();
        Kern_Sigmoid (const Config &c);
        std::string label = "Kern_Sigmoid";
        double operator() (const aed_type2& b, const aed_type2& af);
        aed_type2 derivative(const aed_type2& b, const aed_type2& af);
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff);
        std::string get_label();

};
#endif
