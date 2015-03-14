/* ======================================================================
   Atomistica - Interatomic potential library and molecular dynamics code
   https://github.com/Atomistica/atomistica

   Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
   See the AUTHORS file in the top-level Atomistica directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */

/*
 * Dictionary object.
 */ 

#ifndef __PTRDICT_H
#define __PTRDICT_H

#include <stdio.h>

#ifdef __cplusplus
#define BOOL   bool
#else
#define BOOL   _Bool
#endif
#define FALSE  0
#define TRUE   1


#define MAX_NAME  130
#define MAX_DESCRIPTION  4096

/* Section kinds that exist:
 * SK_SECTION: A section which is a must be
 * SK_MODULE:  A module which *can* be activated
 * SK_1TON:    There can be more than one module of the same name
 * Note: This does only have influence on what is written to
 * a file using ptrdict_write. A section is always written
 * a module only if present in the input file.
 */
#define SK_SECTION  0
#define SK_MODULE   1
#define SK_1TON     2


/* Property kinds */
#define PK_INT  0
#define PK_DOUBLE   1
#define PK_BOOL  2
#define PK_STRING 3
#define PK_FORTRAN_STRING 4
#define PK_POINT 5
#define PK_INTPOINT 6
#define PK_ENUM 7
#define PK_ARRAY1D 8
#define PK_ARRAY2D 9
#define PK_ARRAY3D 10
#define PK_LIST 11
#define PK_STRING_LIST 12
#define PK_FORTRAN_STRING_LIST 13
#define PK_INT_LIST 14


typedef struct __property_t {
  int kind;                             /* Kind of property */
  char name[MAX_NAME+1];                /* Name of this property */
  char description[MAX_DESCRIPTION+1];  /* Help string */
  void *ptr;                            /* Pointer to its memory location */

  int tag;                              /* Additional information, for a string its maximum length */
  int tag2;
  int tag3;
  char *tag4;
  int *tag5;

  /* Access information */
  BOOL provided;                        /* Has this property been provided? */

  /* Traversal information */
  struct __section_t *parent;           /* Parent section */
  struct __property_t *next;            /* Next in this list of properties */
} property_t;


typedef void *(*callback_t)(void *);


typedef struct __section_t {
  int kind;                             /* Kind of section */
  char name[MAX_NAME+1];                /* Name of this section */
  int nalias;                           /* Number of alternative names */
  char alias[MAX_NAME+1];               /* Alternative name */
  char description[MAX_DESCRIPTION+1];  /* Help string */

  /* Callback */
  callback_t callback;                  /* Callback upon access */
  void *tag;                            /* Tag for the callback function */
  void *tag2;                           /* Second tag for the callback function */

  /* Access information */
  BOOL provided;                        /* Was this section present in the input file? */
  BOOL *provided_notification;          /* Notify this variable if the module is provided or not. */

  /* Properties */
  property_t *first_property;           /* List of properties */

  /* Traversal information */
  struct __section_t *parent;           /* Parent section */
  struct __section_t *next;             /* Next in this list of sections */
  struct __section_t *first_child;      /* List of children */
} section_t;


typedef struct {
  FILE *f;
  int column, row;                      /* Column and row information while parsing. */
} parser_t;


/*
 * Registration of sections/parameters
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Create a new group */
section_t *ptrdict_register_group(section_t *self, int kind, char *name,
				 char *description, char *alias);

/* Create a new section */
section_t *ptrdict_register_section(section_t *self, char *name,
				   char *description);

/* Create a new module */
section_t *ptrdict_register_module(section_t *self, BOOL *notification,
				  char *name, char *description);

/* Create a new properties */
void ptrdict_register_integer_property(section_t *self, int *ptr, char *name,
				      char *description);
void ptrdict_register_real_property(section_t *self, double *ptr, char *name,
				   char *description);
void ptrdict_register_boolean_property(section_t *self, BOOL *ptr, char *name,
				      char *description);
void ptrdict_register_string_property(section_t *self, char *ptr, int maxlen,
				     char *name, char *description);

void ptrdict_register_point_property(section_t *self, double *ptr, char *name,
				    char *description);
void ptrdict_register_intpoint_property(section_t *self, int *ptr, char *name,
				       char *description);

void ptrdict_register_enum_property(section_t *self, int *ptr, int nchoices,
				   int lenchoice, char *choices, char *name,
				   char *description);

void ptrdict_register_list_property(section_t *self, double *ptr, int maxlen,
				   int *len, char *name, char *description);
void ptrdict_register_string_list_property(section_t *self, char *ptr,
					  int strlen, int maxlen, int *len,
					  char *name, char *description);
void ptrdict_register_integer_list_property(section_t *self, double *ptr,
					   int maxlen, int *len, char *name,
					   char *description);

void ptrdict_register_array2d_property(section_t *self, double *ptr, int nx,
				      int ny, char *name, char *description);
void ptrdict_register_array3d_property(section_t *self, double *ptr, int nx,
				      int ny, int nz, char *name,
				      char *description);

/* Clean up, remove everything from memory */
void ptrdict_cleanup(section_t *root);



/*
 * Input/output
 */

/* Read ptrdicturation from a stream. */
void ptrdict_from_stream(section_t *root, FILE *f);

/* Read ptrdicturation from a file. */
void ptrdict_read(section_t *root, char *fn);

#ifndef __APPLE__
/* Read ptrdicturation from a string. */
void ptrdict_from_string(section_t *root, char *s);
#endif

/* Write current ptrdicturation to a file. */
void ptrdict_write(section_t *root, char *fn);

#ifdef __cplusplus
};
#endif

#endif
