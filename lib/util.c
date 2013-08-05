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

/*
 * Collection of miscellaneous functions shared throughout source.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "util.h"

#include <stdio.h>


void 
xfree(void *p)
{
    if (p != NULL) {
        free(p);
    }
}

void *
xmalloc(size_t size)
{
    register void *value = malloc(size);
    if (value == NULL) {
        perror("virtual memory exhausted");
        abort();
    }
    return value;
}

void *
xcalloc(size_t count, size_t eltsize)
{
    register void *value = calloc(count, eltsize);
    if (value == NULL) {
        perror("virtual memory exhausted");
        abort();
    }
    return value;
}

void *
xrealloc(void *ptr, size_t size)
{
    register void *value = realloc(ptr, size);
    if (value == NULL) {
        perror("Virtual memory exhausted");
        abort();
    }
    return value;
}



