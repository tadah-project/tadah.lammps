#include "catch2/catch.hpp"
#include "../../CORE/config/config.h"
#include "../../MODELS/functions/f_all.h"
#include "../../MD/m_blr_md.h"


TEST_CASE( "Testing M_BLR_MD Constructor 1", "[M_BLR_MD(bf, config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    BF_Linear bf(c1);
    M_BLR_MD<> m(bf, c1);
}
TEST_CASE( "Testing M_BLR_MD Constructor 2", "[M_BLR_MD(bf, config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    Kern_Linear bf(c1);
    REQUIRE_THROWS(M_BLR_MD<>(bf, c1));
}

TEST_CASE( "Testing M_BLR_MD Constructor 3", "[M_BLR_MD(config)]" ) {
    Config c1("tests_data/valid_configs/valid_config5");
    M_BLR_MD<BF_Linear> m(c1);
}
