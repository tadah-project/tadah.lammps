#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include <string>
#include <map>
#include "element.h"

struct PeriodicTable {
    public:
        PeriodicTable() {
            fill();
        };
        const Element &find_by_name(const std::string &s) const {
            return name.at(s);
        }
        const Element &find_by_symbol(const char *s) const {
            return symbol.at(s);
        }
        const Element &find_by_Z(const int &s) const {
            return Z.at(s);
        }

        void add(const Element e) {
            symbol[e.symbol] = e;
            Z[e.Z] = e;
        }
        void add(const char *symbol_, const char *name_, const int Z_) {
            Element e(symbol_, Z_);
            symbol[symbol_] = e;
            name[name_] = e;
            Z[Z_] = e;
        }
    private:
  struct cmp_str
  {
    bool operator()(char const *a, char const *b) const
    {
      return a[0]==b[0] && a[1]==b[1];
    }
  };
        std::map<const char*, Element, cmp_str> symbol;   // symbol:Element
        std::map<std::string, Element> name;   // name:Element
        std::map<int, Element> Z;   // Z:Element
        void fill() {
            add("H","Hydrogen",1);
            add("He","Helium",2);
            add("Li","Lithium",3);
            add("Be","Beryllium",4);
            add("B","Boron",5);
            add("C","Carbon",6);
            add("N","Nitrogen",7);
            add("O","Oxygen",8);
            add("F","Fluorine",9);
            add("Ne","Neon",10);
            add("Na","Sodium",11);
            add("Mg","Magnesium",12);
            add("Al","Aluminium",13);
            add("Si","Silicon",14);
            add("P","Phosphorus",15);
            add("S","Sulfur",16);
            add("Cl","Chlorine",17);
            add("Ar","Argon",18);
            add("K","Potassium",19);
            add("Ca","Calcium",20);
            add("Sc","Scandium",21);
            add("Ti","Titanium",22);
            add("V","Vanadium",23);
            add("Cr","Chromium",24);
            add("Mn","Manganese",25);
            add("Fe","Iron",26);
            add("Co","Cobalt",27);
            add("Ni","Nickel",28);
            add("Cu","Copper",29);
            add("Zn","Zinc",30);
            add("Ga","Gallium",31);
            add("Ge","Germanium",32);
            add("As","Arsenic",33);
            add("Se","Selenium",34);
            add("Br","Bromine",35);
            add("Kr","Krypton",36);
            add("Rb","Rubidium",37);
            add("Sr","Strontium",38);
            add("Y","Yttrium",39);
            add("Zr","Zirconium",40);
            add("Nb","Niobium",41);
            add("Mo","Molybdenum",42);
            add("Tc","Technetium",43);
            add("Ru","Ruthenium",44);
            add("Rh","Rhodium",45);
            add("Pd","Palladium",46);
            add("Ag","Silver",47);
            add("Cd","Cadmium",48);
            add("In","Indium",49);
            add("Sn","Tin",50);
            add("Sb","Antimony",51);
            add("Te","Tellurium",52);
            add("I","Iodine",53);
            add("Xe","Xenon",54);
            add("Cs","Caesium",55);
            add("Ba","Barium",56);
            add("La","Lanthanum",57);
            add("Ce","Cerium",58);
            add("Pr","Praseodymium",59);
            add("Nd","Neodymium",60);
            add("Pm","Promethium",61);
            add("Sm","Samarium",62);
            add("Eu","Europium",63);
            add("Gd","Gadolinium",64);
            add("Tb","Terbium",65);
            add("Dy","Dysprosium",66);
            add("Ho","Holmium",67);
            add("Er","Erbium",68);
            add("Tm","Thulium",69);
            add("Yb","Ytterbium",70);
            add("Lu","Lutetium",71);
            add("Hf","Hafnium",72);
            add("Ta","Tantalum",73);
            add("W","Tungsten",74);
            add("Re","Rhenium",75);
            add("Os","Osmium",76);
            add("Ir","Iridium",77);
            add("Pt","Platinum",78);
            add("Au","Gold",79);
            add("Hg","Mercury",80);
            add("Tl","Thallium",81);
            add("Pb","Lead",82);
            add("Bi","Bismuth",83);
            add("Po","Polonium",84);
            add("At","Astatine",85);
            add("Rn","Radon",86);
            add("Fr","Francium",87);
            add("Ra","Radium",88);
            add("Ac","Actinium",89);
            add("Th","Thorium",90);
            add("Pa","Protactinium",91);
            add("U","Uranium",92);
            add("Np","Neptunium",93);
            add("Pu","Plutonium",94);
            add("Am","Americium",95);
            add("Cm","Curium",96);
            add("Bk","Berkelium",97);
            add("Cf","Californium",98);
            add("Es","Einsteinium",99);
            add("Fm","Fermium",100);
            add("Md","Mendelevium",101);
            add("No","Nobelium",102);
            add("Lr","Lawrencium",103);
            add("Rf","Rutherfordium",104);
            add("Db","Dubnium",105);
            add("Sg","Seaborgium",106);
            add("Bh","Bohrium",107);
            add("Hs","Hassium",108);
            add("Mt","Meitnerium",109);
            add("Ds","Darmstadtium",110);
            add("Rg","Roentgenium",111);
            add("Cn","Copernicium",112);
            add("Uut","Ununtrium",113);
            add("Fl","Flerovium",114);
            add("Uup","Ununpentium",115);
            add("Lv","Livermorium",116);
            add("Uus","Ununseptium",117);
            add("Uuo","Ununoctium",118);

        }
};
#endif
