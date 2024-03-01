#ifndef REGISTRY_H
#define REGISTRY_H
#include <map>
#include <functional>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>

template< typename T, typename ... Args>
using Factory = std::function< T* (Args ...) >;

template< typename Base, typename ... Args >
struct Registry {

    using Map = std::map<std::string,Factory<Base,Args...>>;
    static Map registry;

    template< typename T >
        struct Register {
            Register( const std::string& name ) {
                registry[ name ] = []( Args... args ) -> T* { return new T(args...); };
            }
        };
};

template< typename Base, typename ... Args >
Base* factory(const std::string &type, Args ... args ) {
    auto it = Registry<Base,Args...>::registry.find( type );
    if ( it!=Registry<Base,Args...>::registry.end() ) {
        return (it->second)(args...);
    }

    std::ostringstream oss;
    oss << "Type: " << type << " of base type: "<< typeid(Base).name()
        << " does not exist." << std::endl;
    throw std::runtime_error(oss.str());
    // return nullptr;
}

#endif
