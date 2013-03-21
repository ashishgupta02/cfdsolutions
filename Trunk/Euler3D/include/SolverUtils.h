/*******************************************************************************
 * File:        SolverUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

using namespace std;

#ifndef _SOLVERUTILS_H
#define	_SOLVERUTILS_H

/*!
 * \class CCreateMap
 * \brief creates a map from a list by overloading operator()
 * \tparam T - the key type in the map
 * \tparam U - the mapped value type in the map
 * \author Ashish Gupta
 *
 * We need this to create static const maps that map strings to enum
 * types.
 */
template <typename T, typename U>
class CCreateMap {
private:
    std::map<T, U> m_map;
public:

    CCreateMap(const T& key, const U& val) {
        m_map[key] = val;
    }

    CCreateMap<T, U>& operator()(const T& key, const U& val) {
        m_map[key] = val;
        return *this;
    }

    operator std::map<T, U>() {
        return m_map;
    }
};

/*!
 * \brief add an enum-based option to the param map and set its default value
 * \param[in] option_name - name of the option as it appears the file
 * \param[in,out] option - the option to associate with option_name
 * \param[in] Tmap - a map from strings to type Tenum
 */
template <class T, class Tenum>
void GetOptionNameValue(const string & option_name, T & option, const map<string, Tenum> & Tmap) {
    typename map<string, Tenum>::const_iterator it;
    it = Tmap.find(option_name);
    if (it == Tmap.end()) {
        cerr << "ERROR: Cannot find " << option_name << " in given map."
                << endl;
        throw (-1);
    }
    option = it->second;
}

/*!
 * \brief add an enum-based option to the param map and set its default value
 * \param[in] option_name - name of the option as it appears the file
 * \param[in,out] option - the option to associate with option_name
 * \param[in] Tmap - a map from strings to type Tenum
 */
template <class T, class Tenum>
void GetOptionValueName(string & option_name, T & option, const map<string, Tenum> & Tmap) {
    typename map<string, Tenum>::const_iterator it;
    option_name.assign("NOT FOUND");
    for (it = Tmap.begin(); it != Tmap.end(); it++)
        if (option == it->second)
            option_name.assign(it->first);
}

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in,out] str - string we want to convert
 */
inline void StringToUpperCase(string & str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in] str - string we want a copy of converted to uppercase
 * \returns a copy of str in uppercase
 */
inline string StringToUpperCase(const string & str) {
    string upp_str(str);
    std::transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::toupper);
    return upp_str;
}

#endif	/* _SOLVERUTILS_H */

