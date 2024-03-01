#ifndef NORMALISER_CORE_H
#define NORMALISER_CORE_H

#include <limits>
#include <stdexcept>
#include "config/config.h"

// We do calculate and store st. dev. and mean
// for the first column of Phi matrix
// for linear kernel, first column are 1s
// so we have to deal with them here
class Normaliser_Core {
    public:
        bool bias;
        int verbose;

        v_type std_dev;
        v_type mean;

        Normaliser_Core () {};
        Normaliser_Core (Config &c):
            bias(c.get<bool>("BIAS")),
            verbose(c.get<int>("VERBOSE"))
    {
        if (c.template get<bool>("NORM")) {
            if (c.exist("NSTDEV")) {
                std_dev.resize(c.size("NSTDEV"));
                mean.resize(c.size("NMEAN"));
                bias = c.template get<bool>("BIAS");
                c.template get<v_type>("NSTDEV",std_dev);
                c.template get<v_type>("NMEAN",mean);
            }
        }

    };

        /** Normalise AED */
        void normalise(aed_type2 &aed) {
            size_t n=aed.size();
            if (bias) {
                for (size_t i=1;i<n;++i) {
                    aed[i]=(aed[i]-mean[i])/std_dev[i];
                }
                //aed.tail(n-1) = (aed.tail(n-1) - mean.tail(n-1)).array() / std_dev.array().tail(n-1);
            }
            else {
                for (size_t i=0;i<n;++i) {
                    aed[i]=(aed[i]-mean[i])/std_dev[i];
                }
                //aed = (aed - mean).array() / std_dev.array();
            }
        }

        /** Normalise FD */
        void normalise(fd_type &fd) {
            // fd is nx3
            size_t n=fd.rows();
            if (bias) {
                for (size_t k=0; k<3; ++k) {
                    for (size_t i=1;i<n;++i) {
                        fd[k](i)=fd(k)(i) / std_dev[i];
                    }
                }
            }
            else {
                for (size_t k=0; k<3; ++k) {
                    for (size_t i=0;i<n;++i) {
                        fd[k](i)=fd(k)(i) / std_dev[i];
                    }
                }
            }
        }
        void normalise(fd_type &fd, size_t k) {
            // fd is nx3
            size_t n=fd.rows();
            if (bias) {
                for (size_t i=1;i<n;++i) {
                    fd[k](i)=fd(k)(i) / std_dev[i];
                }
            }
            else {
                for (size_t i=0;i<n;++i) {
                    fd[k](i)=fd(k)(i) / std_dev[i];
                }
            }
        }

};
#endif
