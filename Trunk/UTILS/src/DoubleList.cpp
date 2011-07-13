/*******************************************************************************
 * File:        DoubleList.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdio.h>
#include "DoubleList.h"

int DoubleList::Redimension(int size) {
    if (size <= 0) {
        if (list != NULL)
            free(list);
        if (data != NULL)
            free(data);
        list = NULL;
        data = NULL;
        dim = max = 0;
        return (1);
    } else if (size != dim) {
        if (list != NULL)
            list = (int*) realloc((void*) list, size * sizeof (int));
        else
            list = (int*) malloc(size * sizeof (int));

        if (data != NULL)
            data = (double*) realloc((void*) data, size * sizeof (double));
        else
            data = (double*) malloc(size * sizeof (double));
        
        if ((list == NULL) || (data == NULL)) {
            fprintf(stderr, "DoubleList:Redimension -- Could not allocate space for list.");
            return (0);
        }
        //for (int i=max; i < size; i++)
        //  list[i] = -1;
        dim = size;
        if (dim < max)
            max = dim;
    }
    return (1);
}

int DoubleList::Is_In_List(int n) {
    int i, y;

    for (y = i = 0; i < max && !y; i++)
        if (list[i] == n)
            y = 1;

    return (y);
}

int DoubleList::Times_In_List(int n) {
    int i, y;

    for (y = i = 0; i < max; i++)
        if (list[i] == n)
            y++;

    return (y);
}

int DoubleList::Index(int n) {
    int i, j;

    j = -1;
    for (i = 0; i < max && j < 0; i++)
        if (list[i] == n)
            j = i;

    return (j);
}

int DoubleList::Check_List(int n, double val) {
    const int INC = 5;
    int new_dim;

    new_dim = dim + INC;
    if (!Is_In_List(n)) {
        if (max >= dim) {
            if (!Redimension(new_dim)) {
                fprintf(stderr, "DoubleList:Check_List -- Could not add to list.");
                return (0);
            }
        }
        list[max]   = n;
        data[max++] = val;
    }
    return (1);
}

int DoubleList::Add_To_List(int n, double val) {
    const int INC = 5;
    int new_dim;

    new_dim = dim + INC;
    if (max < dim) {
        list[max]   = n;
        data[max++] = val;
    } else if (Redimension(new_dim)) {
        list[max]   = n;
        data[max++] = val;
        return (1);
    } else {
        fprintf(stderr,"DoubleList:Add_To_List -- Could not add to list.");
        return (0);
    }
    return (1);
}

int DoubleList::Delete_From_List(int n) {
    int i, j, flag = 0;

    for (i = 0; i < max; i++)
        if (list[i] == n) {
            for (j = i; j < max - 1; j++) {
                list[j] = list[j + 1];
                data[j] = data[j + 1];
            }
            max--;
            flag = 1;
            break;
        }

    return (flag);
}

int DoubleList::Replace(int m, int n, double val) {
    int i, flag = 0;

    for (i = 0; i < max; i++)
        if (list[i] == m) {
            list[i] = n;
            data[i] = val;
            flag = 1;
        }

    return (flag);
}

int DoubleList::Reset(int size) {
    if (size < 0) {
        if (list != NULL)
            free(list);
        if (data != NULL)
            free(data);
        list = NULL;
        data = NULL;
        dim = max = 0;
        return (1);
    } else if (size <= dim) {
        max = 0;
        return (1);
    } else {
        if (list != NULL)
            free(list);
        if (data != NULL)
            free(data);
        list = (int*) malloc(size * sizeof (int));
        data = (double*) malloc(size * sizeof (double));
        if ((list == NULL) || (data == NULL)) {
            fprintf(stderr, "DoubleList:Reset -- Could not allocate space for list.");
            dim = max = 0;
            return (0);
        }
        dim = size;
        max = 0;
        //for (int i=max; i < size; i++)
        //  list[i] = -1;
    }
    return (1);
}

int DoubleList::RemoveDuplicates() {
    if (max > 0) {
        int i, j;
        int value;
        int *map = NULL;
        map = (int*) malloc(max*sizeof(int));
        if (map == NULL) {
            fprintf(stderr, "DoubleList:RemoveDuplicates -- Could not allocate space for list.");
            dim = max = 0;
            return (0);
        }
        for (i = 0; i < max; i++)
            map[i] = 0;
        for (i = 0; i < max; i++) {
            value = list[i];
            for (j = i+1; j < max; j++) {
                if (value == list[j])
                    map[j] = -1;
            }
        }
        value = 1;
        for (i = 1; i < max; i++) {
            if (!map[i]) {
                list[value]   = list[i];
                data[value++] = data[i];
            }
        }
        if (max > value) {
            printf("DoubleList:RemoveDuplicates -- Check for ERROR\n");
            max = value;
        }
        free(map);
    }
    return(1);
}

