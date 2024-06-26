#ifndef STORAGE_H
#define STORAGE_H
#include <stdlib.h>
#include "sph.h"

struct particle {
	int id;
	double radius;
	double volume;
	struct s basic;
	double a[3];
};

struct cell {
	int id;
	double length;
	double center[3];
	int pcount;
	size_t pmemory;
	struct particle **particles;
};

extern int pid;
extern int pcount;
extern size_t pmemory;
extern struct particle *particles;

extern int cid;
extern double clength;
extern int ccount[3];
extern size_t cmemory[3];
extern struct cell ***cells;

extern size_t smemory;

#endif
