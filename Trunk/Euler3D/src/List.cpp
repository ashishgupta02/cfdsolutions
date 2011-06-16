/*******************************************************************************
 * File:        List.cpp
 * Author:      Ashish Gupta
 * Revision:    2
 ******************************************************************************/

#include <stdio.h>
#include "List.h"

int List::Redimension(int size) {
    if (size <= 0) {
        if (list != NULL)
            free(list);
        list = NULL;
        dim = max = 0;
        return (1);
    } else if (size != dim) {
        if (list != NULL)
            list = (int*) realloc((void*) list, size * sizeof (int));
        else
            list = (int*) malloc(size * sizeof (int));
        if (list == NULL) {
            fprintf(stderr, "List:Redimension -- Could not allocate space for list.");
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

int List::Is_In_List(int n) {
    int i, y;

    for (y = i = 0; i < max && !y; i++)
        if (list[i] == n)
            y = 1;

    return (y);
}

int List::Times_In_List(int n) {
    int i, y;

    for (y = i = 0; i < max; i++)
        if (list[i] == n)
            y++;

    return (y);
}

int List::Index(int n) {
    int i, j;

    j = -1;
    for (i = 0; i < max && j < 0; i++)
        if (list[i] == n)
            j = i;

    return (j);
}

int List::Check_List(int n) {
    const int INC = 5;
    int new_dim;

    new_dim = dim + INC;
    if (!Is_In_List(n)) {
        if (max >= dim) {
            if (!Redimension(new_dim)) {
                fprintf(stderr, "List:Check_List -- Could not add to list.");
                return (0);
            }
        }
        list[max++] = n;
    }
    return (1);
}

int List::Add_To_List(int n) {
    const int INC = 5;
    int new_dim;

    new_dim = dim + INC;
    if (max < dim)
        list[max++] = n;
    else if (Redimension(new_dim)) {
        list[max++] = n;
        return (1);
    } else {
        fprintf(stderr,"List:Add_To_List -- Could not add to list.");
        return (0);
    }
    return (1);
}

int List::Delete_From_List(int n) {
    int i, j, flag = 0;

    for (i = 0; i < max; i++)
        if (list[i] == n) {
            for (j = i; j < max - 1; j++)
                list[j] = list[j + 1];
            max--;
            flag = 1;
            break;
        }

    return (flag);
}

int List::Replace(int m, int n) {
    int i, flag = 0;

    for (i = 0; i < max; i++)
        if (list[i] == m) {
            list[i] = n;
            flag = 1;
        }

    return (flag);
}

int List::Reset(int size) {
    if (size < 0) {
        if (list != NULL)
            free(list);
        list = NULL;
        dim = max = 0;
        return (1);
    } else if (size <= dim) {
        max = 0;
        return (1);
    } else {
        if (list != NULL)
            free(list);
        list = (int*) malloc(size * sizeof (int));
        if (list == NULL) {
            fprintf(stderr, "List:Reset -- Could not allocate space for list.");
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

int List::RemoveDuplicates() {
    if (max > 0) {
        int i, j;
        int value;
        int *map = NULL;
        map = (int*) malloc(max*sizeof(int));
        if (map == NULL) {
            fprintf(stderr, "List:RemoveDuplicates -- Could not allocate space for list.");
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
            if (!map[i])
                list[value++] = list[i];
        }
        if (max > value) {
            printf("List:RemoveDuplicates -- Check for ERROR\n");
            max = value;
        }
        free(map);
    }
    return(1);
}


