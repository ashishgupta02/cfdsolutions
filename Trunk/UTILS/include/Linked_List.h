/*******************************************************************************
 * File:        Linked_List.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#ifndef _LINKED_LIST_H
#define	_LINKED_LIST_H

class Linked_Node {
public:
    int data;
    Linked_Node *next;
public:
    Linked_Node(int i, Linked_Node *n) {
        data = i, next = n;
    }
};

class Linked_List {
public:
    Linked_Node *head;
public:
    Linked_List();
    ~Linked_List();
    
    void Insert(int i);
    void Remove(int i);
    int Length();
    int Replace(int i, int j);
    int In_list(int i);
};

#endif	/* _LINKED_LIST_H */

