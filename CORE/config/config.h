#ifndef CONFIG_H
#define CONFIG_H

#include "../utils/utils.h"	// for vector printing
#include "../lib/toml11/toml.hpp"

#include <type_traits>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <limits>


/** \brief A dictionary for key - value pairs.
 *
 * This dictionary object must be used in the workflow
 * to configure other objects in the training process.
 *
 * A Config object stores information such as:
 *  - cutoff distances
 *  - descriptor parameters
 *  - paths to training database files
 *  - and many more...
 *
 * Every object using a Config class specifies a list
 * of required and optional KEY - VALUE pairs.
 *
 * The config keys and values are case insensitive.
 *
 * See documentation for the list of supported KEYS.
 *
 * Below are so called internal keys which are not supposed to be set
 * by the user.
 *
 * <hr>
 * <h1>List of internally used keys:</h1>
 *
 * \warning
 *  \anchor INTERNAL_KEY
 *  **INTERNAL_KEY**\n
 *  Do not try to set internal keys.
 *  Corresponding values are calculated by the library when required.
 *  User can access internal key-value pairs as usually.
 *
 * \anchor DSIZE
 * - <\ref INTERNAL_KEY> **DSIZE** <int>
 *
 *  Total dimension of a descriptor.
 *  DSIZE = SIZE2B + SIZE3B + SIZEMB.
 *  Set by a Calculator.
 *
 * \anchor SIZE2B
 * - <\ref INTERNAL_KEY> **SIZE2B** <int>
 *
 *  Dimension of a two-body descriptor.
 *  Set by a Calculator.
 *
 * \anchor SIZE3B
 * - <\ref INTERNAL_KEY> **SIZE3B** <int>
 *
 *  Dimension of a three-body descriptor.
 *  Set by a Calculator.
 *
 * \anchor SIZEMB
 * - <\ref INTERNAL_KEY> **SIZEMB** <int>
 *
 *  Dimension of a many-body descriptor.
 *  Set by a Calculator.
 *
 * \anchor ATOMS
 * - **ATOMS** <string> A1 A2 ...
 *
 *  List of atom types. Read from training database files.
 *  e.g. `ATOMS Ti Nb`
 *
 * \anchor RCUTMAX
 * - **RCUTMAX** <string> A1 A2 ...
 *
 *  Maximum cutoff distance, max(RCUT2B,RCUT3B,RCUTMB)
 *
 * \anchor ESTDEV
 * - **ESTDEV** <double> N
 *
 * 	Energy standard deviation (energy/atom).
 *
 * \anchor FSTDEV
 * - **FSTDEV** <double> N
 *
 * 	Force standard deviation (energy/distance).
 *
 * \anchor SSTDEV
 * - **SSTDEV** <double> xx xy xz yy yz zz
 *
 * 	Stress standard deviations (energy units).
 *
 */
class Config {

    public:
        /** \brief Create object with default values. */
        Config();

        /** \brief Create object and read config values from the file.
         *
         * The values from the file take precedence over config defaults.
         */
        Config (std::string fn);

        /** \brief Return all stored values for a key as a vector
         *
         * Usage example:
         *
         * \code{.cpp}
         * # Create new config object and read config file
         * Config config("config_file");
         *
         * # store value(s) for key "DBFILE" in a vector v
         * auto v = config("DBFILE");
         * \endcode
         *
         */
        const std::vector<std::string>& operator()(const std::string key) const;

        /** \brief Fill array with the values from the key.
         *
         *  The array must be initialised to appropriate size.
         *
         * Usage example:
         *
         * \code{.cpp}
         * # Create new config object and read config file
         * Config config("config_file");
         *
         * # Fill array with config values of type double
         * std::vector<double> vd(10);
         * c.get<std::vector<double>>("RCUT2B", vd);
         * \endcode
         *
         * \tparam T indexable array
         */
        template <typename T>
            void get(const std::string key, T & value) const {
                if (c.count(to_upper(key))==0)
                    throw std::runtime_error(knf+key);
                int i=0;
                for (auto v : c.at(to_upper(key))) {
                    std::istringstream iss(v);
                        iss >> value[i++];
                }
            }

        /** \brief Return the first value for the key.
         *
         * Usage example:
         *
         * \code{.cpp}
         * # Create new config object and read config_file
         * Config config("config_file");
         *
         * # get value from the config file and assign to d
         * double d = c.get<double>("RCUT2B");
         * \endcode
         *
         * \tparam T must be built-in type, e.g. int, double, char
         */
        template <typename T>
            T get(const std::string key) const {
                if (c.count(to_upper(key))==0)
                    throw std::runtime_error(knf+key);
                T value;
                std::istringstream iss(c.at(to_upper(key))[0]);
                iss >> std::boolalpha >> value;
                return value;
            }

        /** \brief Return n-th value for a key with a specified type
         *
         * Usage example:
         *
         * \code{.cpp}
         * # Create new config object and read config_file
         * Config config("config_file");
         *
         * # get value from the config file and assign to d
         * double d = c.get<double>("CGRID2B",3);
         * \endcode
         *
         * \tparam T must be built-in type, e.g. int, double, char
         */
        template <typename T>
            T get(const std::string key,size_t n) const {
                if (c.count(key)==0)
                    throw std::runtime_error(knf+key);
                if (c.at(key).size()-1<n)
                    throw std::runtime_error("Requested value is larger than available at key: "+key);
                T value;
                std::istringstream iss(c.at(key)[n]);
                iss >> std::boolalpha >> value;
                return value;
            }

        /** \brief Add config value to a key.
         *
         * If the key does not exist then new entry is added.
         * If the key exists and allow for multiple values then the value is appended.
         * If the key exists and and is full then throw.
         *
         * This method is case insensitive, i.e.,
         * added keys are always converted to upper case.
         *
         * The value is converted to string.
         *
         * This method allow to add also internal keys.
         */
        template <typename T>
            void add(const std::string key, T val) {
                int n=get<int>("OUTPREC");
                std::string k = to_upper(key);
                if (!tdc.contains(k) && internal_keys.count(k)==0) {
                    throw std::runtime_error(
                            "This key is not allowed: "+k+"\n");
                }
                const auto& allowed_key = toml::find(tdc, k);
                int valN  = toml::find<int>(allowed_key, "value_N",0);
                if (internal_keys.count(k)>0) {
                    if (c[k].size() >= (size_t)abs(internal_keys.at(k))) {
                        throw std::runtime_error(
                                "Internal Config key is already set/full: "+k+"\n");
                    }
                }
                else if (c[k].size() >= (size_t)abs(valN)) {
                    throw std::runtime_error("Config key is full, KEY: "
                            +k+" VALUE: "+to_string(val,n)+"\n");
                }
                std::string val_str = to_string(val,n);
                std::istringstream iss(val_str);
                std::string value;
                while (iss >> std::boolalpha >> value) {
                    c[k].push_back(to_string(value,n));
                }
            }

        /** Remove the key from the config.
         *
         * Return number of keys removed, i.e., 0 or 1
         */
        size_t remove(const std::string key);

        /** \brief Add key-value pairs to the existing config object from the file.
         *
         * const std::string fn is a filename containing config.
         */
        void add(const std::string fn);

        /** \brief Return a number of values stored by the key.
         *
         * Throws if the key does not exists.
         */
        size_t size(const std::string key) const;

        /** \brief Return true if the key is in this object. */
        bool exist(const std::string key) const;

        void postprocess();

        /** \brief Check is this object configured for training.
         *
         * This helper method allows to verify whether the config is
         * suitable for trainig. It checks some basic properties
         * but does not cover every possibility.
         *
         * TODO return false instead of throwing?
         */
        bool check_for_training();

        /** \brief Check is this object configured for prediction.
         *
         * TODO Not implemented
         */
        bool check_for_predict();

        /** \brief Verify is this object configured as a potential file.
         *
         * TODO Not implemented
         */
        bool check_pot_file();
        void set_defaults();
        void clear_internal_keys();

        /** \brief Return true if both configs are the same. */
        bool operator==(const Config&) const;
    private:
        //toml::table tdc;
        toml::value tdc;
        std::map<std::string, std::vector<std::string>> c;
        void read(std::string fn);

        /** \brief Send all key-value pairs to the stream. */
        friend std::ostream& operator<<(std::ostream& os, const Config& config)
        {
            for (auto pair : config.c) {
                const auto& allowed_key = toml::find(config.tdc, pair.first);
                int valN  = toml::find<int>(allowed_key, "value_N",0);
                if (valN > 1) {
                    for (auto v : pair.second) {
                        os << pair.first << "     " << v << std::endl;
                    }
                }
                else if (valN < 1) {
                    os << pair.first;
                    for (auto v : pair.second) {
                        os << "     " << v;
                    }
                    os << std::endl;
                }
                else {
                    os << pair.first;
                    for (auto v : pair.second) {
                        os << "      " << v;
                    }
                    os << std::endl;
                }
            }
            return os;
        }

        static std::map<std::string,int> create_internal_keys() {
            std::map<std::string,int> m;

            // number indicates max number of allowed values
            // specified in a config file
            // INTERNALLY USED KEYS
            m["DSIZE"]=1;
            m["SIZE2B"]=1;
            m["SIZE3B"]=1;
            m["SIZEMB"]=1;
            m["ATOMS"]=-std::numeric_limits<int>::max();   // TODO Do we need it?
            m["RCUTMAX"]=1;
            m["ESTDEV"]=1;
            m["FSTDEV"]=1;
            m["SSTDEV"]=-6;

            return m;
        }

		int read_tadah_configuration();

        /** Config file allowed keys and their multiplicity
         *
         * abs(N)=1 max one key provided in a config file
         * abs(N)>1 more than one key can be provided in a config file:
         *
         * Negative key value indicates how key - value(s) are printed.
         *
         * Positive:
         * DBFILE 1
         * DBFILE 2
         *
         * Negative:
         * WEIGHTS 1 2 3 4
         */
        static const std::map<std::string,int> internal_keys;
        const static std::string knf;//="Key not found: ";
};
#endif
