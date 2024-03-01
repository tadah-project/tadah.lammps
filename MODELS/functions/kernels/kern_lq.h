#ifndef KERN_LQ_H
#define KERN_LQ_H

#include "kern_base.h"
/**
 * Linear + Quadratic kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T \mathbf{y} +
 * \Big(\mathbf{x}^T \mathbf{y} \Big)^2
 * = \sum_i x^{(i)} y^{(i)} + \Big(\sum_i x^{(i)} y^{(i)} \Big)^2
 * \f]
 *
 *  @see Kern_Base Kern_Linear Kern_Quadratic
 */
class Kern_LQ : public virtual Kern_Base {
    public:
        Kern_LQ ();
        Kern_LQ (const Config &c);
        /**
         * Label used for this class
         */
        std::string label = "Kern_LQ";
        double operator() (const aed_type2& b, const aed_type2& af);
        aed_type2 derivative(const aed_type2& b, const aed_type2& af);
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff);
        std::string get_label();

};
#endif
