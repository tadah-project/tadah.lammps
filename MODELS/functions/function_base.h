#ifndef FUNCTION_BASE_h
#define FUNCTION_BASE_h

#include "../../CORE/core_types.h"
#include "../../CORE/config/config.h"


/** Base class for Kernels and Basis Functions */
struct Function_Base {

    private:
        int verbose;

    public:
        Function_Base() {
            verbose=0;
        }
        Function_Base(const Config &c) {
            verbose = c.get<int>("VERBOSE");
        }
        void set_verbose(int v) { verbose=v; }
        int get_verbose() { return verbose; }

        // Derived classes must implement Derived() and Derived(Config)
        virtual ~Function_Base() {};

        // For Kernels
        virtual double operator() (const aed_type2& , const aed_type2& ) { return 0.0; };
        virtual aed_type2 derivative(const aed_type2& , const aed_type2& ) { return aed_type2(); };
        virtual double prime(const aed_type2& , const aed_type2& ,
                const aed_type2& ) { return 0.0; }
        virtual void set_basis(const Matrix ) {};

        /** \brief Return label of this object */
        virtual std::string get_label()=0;
        virtual double epredict(const t_type &, const aed_type2& )=0;
        virtual double fpredict(const t_type &, const fd_type &,
                const aed_type2&, const size_t )=0;
        virtual force_type fpredict(const t_type &, const fd_type &,
                const aed_type2& )=0;
};
//template<> inline Registry<Function_Base>::Map Registry<Function_Base>::registry{};
//template<> inline Registry<Function_Base,Config&>::Map Registry<Function_Base,Config&>::registry{};
#endif
