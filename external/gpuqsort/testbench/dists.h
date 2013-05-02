#ifndef DIST_H
#define DIST_H

typedef unsigned int element;

unsigned long genrand_int32(void);
const char* getDistString(unsigned int dist);
unsigned int* readModel(int& size, char* name);
void dist(element* data, unsigned int size, int type);

#endif
