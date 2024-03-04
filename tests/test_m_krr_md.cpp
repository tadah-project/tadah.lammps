#include "catch2/catch.hpp"
#include "../../CORE/config/config.h"
#include "../../MODELS/functions/f_all.h"
#include "../../MD/m_krr_md.h"


TEST_CASE( "Testing M_KRR_MD Constructor 1", "[M_BLR_MD(kern, config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    BF_Linear kern(c1);
    REQUIRE_THROWS(M_KRR_MD<>(kern, c1));
}
TEST_CASE( "Testing M_KRR_MD Constructor 2", "[M_BLR_MD(kern, config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    Kern_Linear kern(c1);
    M_KRR_MD<> m(kern, c1);
}

TEST_CASE( "Testing M_KRR_MD Constructor 3", "[M_BLR_MD(config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    M_KRR_MD<Kern_Linear> m(c1);
}
