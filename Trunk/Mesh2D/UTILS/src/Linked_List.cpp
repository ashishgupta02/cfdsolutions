/*
 * File:   Linked_List.cpp
 * Author: Ashish Gupta
 *
 * Created on April 20, 2010, 6:45 PM
 */

#include "Linked_List.h"

// *****************************************************************************
// *****************************************************************************
Linked_List::Linked_List() {
    head = 0;
}

// *****************************************************************************
// *****************************************************************************
Linked_List::~Linked_List() {
    Linked_Node *current = head;
    while (current) {
        Linked_Node *nxt = current->next;
        delete current;
        current = nxt;
    }
    head = 0;
}

// *****************************************************************************
// *****************************************************************************
void Linked_List::Insert(int i) {
    head = new Linked_Node(i, head);
}

// *****************************************************************************
// *****************************************************************************
int Linked_List::Length() {
    int l = 0;
    Linked_Node *current = head;
    while (current) {
        l++;
        current = current->next;
    }
    return (l);
}

// *****************************************************************************
// *****************************************************************************
void Linked_List::Remove(int i) {
    Linked_Node *current = head;
    Linked_Node *prev = 0;
    while (current && current->data != i) {
        prev = current;
        current = current->next;
    }
    if (current && current->data == i) {
        Linked_Node *nxt = current->next;
        if (prev == 0)
            head = nxt;
        else
            prev->next = nxt;
        delete current;
    }
}

// *****************************************************************************
// *****************************************************************************
int Linked_List::Replace(int i, int j) {
    int flag = 0;
    Linked_Node *current = head;
    while (current && current->data != i) {
        current = current->next;
    }
    if (current && current->data == i) {
        current->data = j;
        flag = 1;
    }

    return (flag);
}

// *****************************************************************************
// *****************************************************************************
int Linked_List::In_list(int i) {
    int flag = 0;
    Linked_Node *current = head;
    while (current && current->data != i) {
        current = current->next;
    }
    if (current && current->data == i)
        flag = 1;

    return (flag);
}

