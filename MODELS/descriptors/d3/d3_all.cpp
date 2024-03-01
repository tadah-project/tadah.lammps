#include "d3_base.h"
#include "d3_dummy.h"

template<> Registry<D3_Base,Config&>::Map Registry<D3_Base,Config&>::registry{};

Registry<D3_Base,Config&>::Register<D3_Dummy> D3_Dummy_1( "D3_Dummy" );
