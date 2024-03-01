#ifndef D_BASE_H
#define D_BASE_H

#include "../../CORE/config/config.h"
#include "../../CORE/core_types.h"
#include "../../CORE/registry.h"

/** \brief Base class for all descriptor types.
 *
 */
class D_Base {
    public:
        static void get_grid(const Config &c, std::string key, v_type &v) {
            if (c.get<double>(key) < 0) {
                // generate grid automatically if first element
                // is<zero and int. If it is double than assume that
                // grid is provided manually and min value is just smaller
                // than zero.
                double int_part;
                if (std::modf ( c.get<double>(key), &int_part ) == 0.0) {
                    // automatic generation
                    size_t cg = abs(c.get<int>(key,1));
                    double cgmin = c.get<double>(key,2);
                    double cgmax = c.get<double>(key,3);
                    if (c.get<int>(key) == -1) {
                        v = linspace(cgmin, cgmax, cg);
                    }
                    else if (c.get<int>(key) == -2) {
                        v = logspace(cgmin, cgmax, cg, exp(1));
                    }
                    else {
                        throw std::runtime_error(key+" algorithm not supported\n");
                    }
                }
                else {
                    v.resize(c.size(key));
                    c.get<v_type>(key,v);
                }
            }
            else {
                v.resize(c.size(key));
                c.get<v_type>(key,v);
            }
        }
        int verbose;
        virtual ~D_Base() {};

        /** \brief Return dimension of the descriptor.
         */
        virtual size_t size()=0;

        /** \brief Return label of this descriptor.
         */
        virtual std::string label()=0;
};
#endif
