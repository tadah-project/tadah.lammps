#include "f_all.h"
#include "../../CORE/registry.h"

template<> Registry<Function_Base>::Map Registry<Function_Base>::registry{};
template<> Registry<Function_Base,Config&>::Map Registry<Function_Base,Config&>::registry{};

Registry<Function_Base>::Register<BF_Linear> BF_Linear_1( "BF_Linear" );
Registry<Function_Base,Config&>::Register<BF_Linear> BF_Linear_2( "BF_Linear" );
Registry<Function_Base>::Register<BF_Polynomial2> BF_Polynomial2_1( "BF_Polynomial2" );
Registry<Function_Base,Config&>::Register<BF_Polynomial2> BF_Polynomial2_2( "BF_Polynomial2" );
Registry<Function_Base>::Register<Kern_Linear> Kern_Linear_1( "Kern_Linear" );
Registry<Function_Base,Config&>::Register<Kern_Linear> Kern_Linear_2( "Kern_Linear" );
Registry<Function_Base>::Register<Kern_LQ> Kern_LQ_1( "Kern_LQ" );
Registry<Function_Base,Config&>::Register<Kern_LQ> Kern_LQ_2( "Kern_LQ" );
Registry<Function_Base>::Register<Kern_Polynomial> Kern_Polynomial_1( "Kern_Polynomial" );
Registry<Function_Base,Config&>::Register<Kern_Polynomial> Kern_Polynomial_2( "Kern_Polynomial" );
Registry<Function_Base>::Register<Kern_Quadratic> Kern_Quadratic_1( "Kern_Quadratic" );
Registry<Function_Base,Config&>::Register<Kern_Quadratic> Kern_Quadratic_2( "Kern_Quadratic" );
Registry<Function_Base>::Register<Kern_RBF> Kern_RBF_1( "Kern_RBF" );
Registry<Function_Base,Config&>::Register<Kern_RBF> Kern_RBF_2( "Kern_RBF" );
Registry<Function_Base>::Register<Kern_Sigmoid> Kern_Sigmoid_1( "Kern_Sigmoid" );
Registry<Function_Base,Config&>::Register<Kern_Sigmoid> Kern_Sigmoid_2( "Kern_Sigmoid" );

