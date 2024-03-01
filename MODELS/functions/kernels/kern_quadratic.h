#ifndef KERN_QUADRATIC_H
#define KERN_QUADRATIC_H
#include "kern_base.h"
/**
 * Quadratic kernel - special case of 2nd order polynomial kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \Big( \mathbf{x}^T \mathbf{y} \Big)^2
 * = \Big( \sum_i x^{(i)} y^{(i)} \Big)^2
 * \f]
 *
 *  @see Kern_Base
 */
class Kern_Quadratic : public virtual Kern_Base {
    public:
        Kern_Quadratic ();
        Kern_Quadratic (const Config &c);
        /**
         * Label used for this class
         */
        std::string label = "Kern_Quadratic";
        double operator() (const aed_type2& b, const aed_type2& af);
        aed_type2 derivative(const aed_type2& b, const aed_type2& af);
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff);
        std::string get_label();

};
#endif
