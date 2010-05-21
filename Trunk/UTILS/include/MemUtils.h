/* 
 * File:   MemUtils.h
 * Author: Ashish Gupta
 *
 * Created on April 21, 2010, 3:30 PM
 */

#include <algorithm>
#include <iostream>
#include <cstddef>

#ifndef _MEMUTILS_H
#define	_MEMUTILS_H

template <typename T>
T *Mem_Relloc(T *array, size_t old_size, size_t new_size) {
    T *temp = new T[new_size];

    delete [] array;

    return std::copy(array, array + old_size, temp);
}

#endif	/* _MEMUTILS_H */

