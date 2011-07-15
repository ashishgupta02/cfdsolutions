/*******************************************************************************
 * File:        DoubleList.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _DOUBLELIST_H
#define	_DOUBLELIST_H

#include <stdio.h>
#include <stdlib.h>

class DoubleList {
private:
    int dim;
public:
    int    *list;
    double *data;
    int    max;

    // Constructor
    DoubleList() {
        dim = max = 0;
        list = NULL;
        data = NULL;
    }

    // Destructor
    ~DoubleList() {
        dim = max = 0;
        if (list != NULL) {
            free(list);
            list = NULL;
        }
        if (data != NULL) {
            free(data);
            data = NULL;
        }
    }

    void construct() {
        dim = max = 0;
        list = NULL;
        data = NULL;
    }

    void destruct() {
        dim = max = 0;
        if (list != NULL) {
            free(list);
            list = NULL;
        }
        if (data != NULL) {
            free(data);
            data = NULL;
        }
    }

    // Copy Constructor
    DoubleList(DoubleList* s) {
        int i;
        Redimension(s->dim);
        max = 0;
        for (i = 0; i < s->max; i++)
            Add_To_List(s->list[i], s->data[i]);
    }

    // Copy constructor
    DoubleList(DoubleList& s) {
        int i;
        Redimension(s.dim);
        max = 0;
        for (i = 0; i < s.max; i++)
            Add_To_List(s.list[i], s.data[i]);
    }

    // operator = for use when LHS already exists
    DoubleList & operator=(DoubleList& s) {
        int i;
        Redimension(s.dim);
        max = 0;
        for (i = 0; i < s.max; i++)
            Add_To_List(s.list[i], s.data[i]);
        return *this;
    }

    void print_dim(FILE *outf) {
        fprintf(outf, "\nInteger-Double list dimension =%d", dim);
    }

    void print(FILE *outf) {
        print_dim(outf);
        fprintf(outf, "\nInteger list maximum index =%d", max);
        for (int i = 0; i < max; i++)
            fprintf(outf, "\nlist(%d)= %d, %lf", i, list[i], data[i]);
    }

    int Redimension(int size);
    int Is_In_List(int n);
    int Check_List(int n, double val);
    int Add_To_List(int n, double val);
    int Delete_From_List(int n);
    int Replace(int m, int n, double val);
    int Times_In_List(int n);
    int Index(int n);
    int Reset(int size);
    int RemoveDuplicates();
};

#endif	/* _DOUBLELIST_H */

