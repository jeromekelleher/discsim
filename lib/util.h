/*
** Copyright (C) 2013 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This file is part of discsim.
** 
** discsim is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** discsim is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with discsim.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>

#define ERR_ALLOC_FAILED -1
#define ERR_BAD_PIXEL_SIZE -2 
#define ERR_OUT_OF_INT_MAP_NODES -3
#define ERR_OUT_OF_SET_MAP_NODES -4 
#define ERR_OUT_OF_AVL_SET_NODES -5 
#define ERR_OUT_OF_INDIVIDUALS -6 
#define ERR_AVL_OP_FAILED -7
#define ERR_MAX_OCCUPANCY_EXCEEDED -8
#define ERR_BAD_DIMENSION -9
#define ERR_EVENT_TOO_LARGE -10

typedef struct {
    double rate;
    double r;
    double u;
} event_class_t;


void xfree(void *p);
void * xmalloc(size_t size);
void * xcalloc(size_t count, size_t eltsize);
void * xrealloc(void *ptr, size_t size);

const char * discsim_error_message(int err);

#endif /*__UTIL_H__*/
