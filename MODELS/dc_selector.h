#ifndef DC_SELECTOR_H
#define DC_SELECTOR_H

#include "../CORE/config/config.h"
#include "descriptors/d_all.h"
#include "cutoffs/cut_all.h"
#include "../CORE/registry.h"

class DC_Selector {
    public:
        Cut_Base *c2b=nullptr;
        Cut_Base *c3b=nullptr;
        Cut_Base *cmb=nullptr;
        D2_Base *d2b=nullptr;
        D3_Base *d3b=nullptr;
        DM_Base *dmb=nullptr;

        Config config;

        DC_Selector ()
        {};

        // Copy constructor
        DC_Selector(const DC_Selector& )
        {
            std::cout << "Copy constructor called\n";
        }
        // Const Assignment operator
        // Instead of deep copying
        // we use a trick where we just copy
        // the config file and run init()
        // as in a costructor
        DC_Selector& operator=(const DC_Selector& dc)
        {
            config=dc.config;
            init();
            return *this;
        }
        // Assignment operator
        DC_Selector& operator=(DC_Selector& dc)
        {
            std::swap(config,dc.config);
            std::swap(c2b,dc.c2b);
            std::swap(c3b,dc.c3b);
            std::swap(cmb,dc.cmb);
            std::swap(d2b,dc.d2b);
            std::swap(d3b,dc.d3b);
            std::swap(dmb,dc.dmb);
            return *this;
        }

        DC_Selector (const Config &c):
            config(c)
    {
        init();
    };
        void init() {
            double rcutzero = 0.0; // for dummies

            size_t bias=0;
            if (config.get<bool>("BIAS"))
                bias++;

            if (config.get<bool>("INIT2B")) {
                double rcut2b = config.get<double>("RCUT2B");
                c2b = factory<Cut_Base,double>( config.get<std::string>("RCTYPE2B"), rcut2b );
                d2b = factory<D2_Base,Config&>( config.get<std::string>("TYPE2B"), config );
                d2b->fidx = bias;
            }
            else {
                c2b = factory<Cut_Base,double>( "Cut_Dummy", rcutzero );
                d2b = factory<D2_Base,Config&>( "D2_Dummy", config );
                d2b->fidx = bias;
            }

            if (config.get<bool>("INIT3B")) {
                double rcut3b = config.get<double>("RCUT3B");
                c3b = factory<Cut_Base,double>( config.get<std::string>("RCTYPE3B"), rcut3b );
                d3b = factory<D3_Base,Config&>( config.get<std::string>("TYPE3B"), config );
                // d3b->fidx = bias;
            }
            else {
                c3b = factory<Cut_Base,double>( "Cut_Dummy", rcutzero );
                d3b = factory<D3_Base,Config&>( "D3_Dummy", config );
                // d3b->fidx = bias;
            }

            if (config.get<bool>("INITMB")) {
                double rcutmb = config.get<double>("RCUTMB");
                cmb = factory<Cut_Base,double>( config.get<std::string>("RCTYPEMB"), rcutmb );
                dmb = factory<DM_Base,Config&>( config.get<std::string>("TYPEMB"), config );
                dmb->fidx = bias + d2b->size();
            }
            else {
                cmb = factory<Cut_Base,double>( "Cut_Dummy", rcutzero );
                dmb = factory<DM_Base,Config&>( "DM_Dummy", config );
                dmb->fidx = bias + d2b->size();
            }
        }
        ~DC_Selector()
        {
            if (d2b)
                delete d2b;
            if (c2b)
                delete c2b;
            if (c3b)
                delete c3b;
            if (d3b)
                delete d3b;
            if (cmb)
                delete cmb;
            if (dmb)
                delete dmb;
        }
};
#endif
