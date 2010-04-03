/* 
 * File:   List.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */


#ifndef _LIST_H
#define	_LIST_H

#include <stdio.h>
#include <stdlib.h>

class List {
private:
    int dim;
public:
    int *list;
    int max;

    // Constructor
    List() {
        dim = max = 0;
        list = NULL;
    }

    // Distructor
    ~List() {
        dim = max = 0;
        if (list != NULL) {
            free(list);
            list = NULL;
        }
    }

    void construct() {
        dim = max = 0;
        list = NULL;
    }

    void destruct() {
        dim = max = 0;
        if (list != NULL) {
            free(list);
            list = NULL;
        }
    }

    // Copy Constructor
    List(List* s) {
        int i;
        Redimension(s->dim);
        max = 0;
        for (i = 0; i < s->max; i++)
            Add_To_List(s->list[i]);
    }

    // Copy constructor
    List(List& s) {
        int i;
        Redimension(s.dim);
        max = 0;
        for (i = 0; i < s.max; i++)
            Add_To_List(s.list[i]);
    }

    // operator = for use when LHS already exists
    List & operator=(List& s) {
        int i;
        Redimension(s.dim);
        max = 0;
        for (i = 0; i < s.max; i++)
            Add_To_List(s.list[i]);
        return *this;
    }

    void print_dim(FILE *outf) {
        fprintf(outf, "\nInteger list dimension =%d", dim);
    }

    void print(FILE *outf) {
        print_dim(outf);
        fprintf(outf, "\nInteger list maximum index =%d", max);
        for (int i = 0; i < max; i++)
            fprintf(outf, "\nlist(%d)= %d", i, list[i]);
    }

    int Redimension(int size);
    int Is_In_List(int n);
    int Check_List(int n);
    int Add_To_List(int n);
    int Delete_From_List(int n);
    int Replace(int m, int n);
    int Times_In_List(int n);
    int Index(int n);
    int Reset(int size);
    int RemoveDuplicates();
};

#endif	/* _LIST_H */

