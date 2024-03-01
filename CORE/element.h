#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <functional>
#include <string>

struct Element {
private:
  void check(const std::string &symbol, const int Z) {
    if (symbol.size() > 3)
      throw std::runtime_error("Wrong symbol: "+symbol);
    if (Z < 1)
      throw std::runtime_error("Atomic number Z < 1");
    if (Z > 118)
      throw std::runtime_error("Atomic number Z > 118");
  }

  friend std::ostream& operator<<(std::ostream& os, const Element& element)
  {
    os << "Symbol: " << element.symbol[0] << element.symbol[1] << std::endl;
    os << "Atomic Number: " << element.Z << std::endl;
    return os;
  }

public:
  char symbol[2];
  int Z;
  Element() {};
  Element(const char *symbol_, const int Z_, bool check_=true)
  {
    if (check_) check(symbol_, Z_);
    symbol[0]=symbol_[0];
    symbol[1]=symbol_[1];
    Z=Z_;
  };
  bool operator==(const Element& other) const {
    return this->symbol == other.symbol && this->Z == other.Z;
  };
  bool operator<(const Element& other) const {
    return this->Z < other.Z;
  };
};

namespace std {
template <>
class hash<Element> {
public:
  size_t operator()(const Element &e) const
  {
    std::hash<int> hash_f;
    return hash_f(e.Z);
  }
};
}
#endif
