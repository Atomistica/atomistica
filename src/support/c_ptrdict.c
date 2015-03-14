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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptrdict.h"

#define MIN(x, y)  (x < y ? x : y)

#define CHECK_NAME(fn, name) \
  if (strlen(name) > MAX_NAME) { \
    printf("["#fn"] Internal error: Section name too long: '%s'.\n", name); \
    exit(1); \
  } while(0)

#define CHECK_DESCRIPTION(fn, name) \
  if (strlen(name) > MAX_DESCRIPTION) { \
    printf("["#fn"] Internal error: Description to long: '%s'.\n", name); \
    exit(1); \
  } while(0)


/* Create a new group */
section_t *ptrdict_register_group(section_t *self, int kind, char *name,
				                          char *description, char *alias)
{
  section_t *new_section;
  section_t *i;

  /* Sanity checks */
  CHECK_NAME(ptrdict_register_section, name);
  CHECK_DESCRIPTION(ptrdict_register_section, description);

#ifdef DEBUG
  printf("ptrdict_register_group: %p %i %s %s\n", self, kind, name, description);
#endif

  if (kind != SK_SECTION && kind != SK_MODULE && kind != SK_1TON) {
    printf("[ptrdict_register_section] Internal error: Unknown section kind: %i.\n", kind);
    exit(1);
  }

  new_section = (section_t*) malloc(sizeof(section_t));
  new_section->kind = kind;
  strcpy(new_section->name, name);
  if (alias) {
    if (alias[0] != 0) {
      strcpy(new_section->alias, alias);
      new_section->nalias = 1;
    } else {
      new_section->nalias = 0;
    }
  } else {
    new_section->nalias = 0;
  }
  strcpy(new_section->description, description);
  new_section->callback = NULL;
  new_section->tag = NULL;
  new_section->tag2 = NULL;
  new_section->provided = FALSE;
  new_section->provided_notification = NULL;

  new_section->first_property = NULL;

  new_section->parent = self;
  new_section->first_child = NULL;
  new_section->next = NULL;

/*   if (self) */
/*     new_section->next = self->first_child; */
/*   else */
/*     new_section->next = NULL; */

/*   /\* If self == NULL this is the root section *\/ */
/*   if (self) */
/*     self->first_child = new_section; */

  /* If self == NULL this is the root section */
  if (self) {

    i = self->first_child;
    if (i) {

      /* Append the new section to the end of the list. */

      while (i->next)  i = i->next;
      i->next = new_section;

    } else {

      /* If the current section does not have a child,
	 this is the first one. */

      self->first_child = new_section;

    };
  };

#ifdef DEBUG
  printf("ptrdict_register_group: new group object @ %p\n", new_section);
#endif

  return new_section;
}


/* Create a new section */
section_t *ptrdict_register_section(section_t *self, char *name,
				                            char *description)
{
  return ptrdict_register_group(self, SK_SECTION, name, description, NULL);
}


/* Create a new module */
section_t *ptrdict_register_module(section_t *self, BOOL *notification,
				                           char *name, char *description)
{
  section_t *s;
  
  s = ptrdict_register_group(self, SK_MODULE, name, description, NULL);
  s->provided_notification = notification;
  *notification = FALSE;

#ifdef DEBUG
  printf("ptrdict_register_module: new module object @ %p\n", s);
#endif

  return s;
}


/* General property */
property_t *ptrdict_register_property(section_t *self, int kind, void *ptr,
				     char *name, char *description)
{
  property_t *new_property, *last;

  /* Sanity checks */
  CHECK_NAME(ptrdict_register_property, name);
  CHECK_DESCRIPTION(ptrdict_register_property, description);

#ifdef DEBUG
  printf("ptrdict_register_property: %p %i %p %s %s\n", self, kind, ptr, name, description);
#endif

  new_property = (property_t*) malloc(sizeof(property_t));
  new_property->kind = kind;
  strcpy(new_property->name, name);
  strcpy(new_property->description, description);
  new_property->ptr = ptr;

  new_property->parent = self;
  new_property->next = NULL;

  new_property->provided = FALSE;

  /* Insert at the end so the list is properly sorted. */
  if (!self->first_property)
    self->first_property = new_property;
  else {
    last = self->first_property;
    while (last->next)
      last = last->next;
    last->next = new_property;
  }

  return new_property;
}


/* Integer property */
void ptrdict_register_integer_property(section_t *self, int *ptr, char *name,
				      char *description)
{
#ifdef DEBUG
  printf("ptrdict_register_integer_property: %p %p\n", self, *self);
#endif

  ptrdict_register_property(self, PK_INT, ptr, name, description);
}


/* Double property */
void ptrdict_register_real_property(section_t *self, double *ptr, char *name,
				   char *description)
{
#ifdef DEBUG
  printf("ptrdict_register_real_property: %p %p\n", self, *self);
#endif

  ptrdict_register_property(self, PK_DOUBLE, ptr, name, description);
}


/* Boolean property */
void ptrdict_register_boolean_property(section_t *self, BOOL *ptr, char *name,
				      char *description)
{
  ptrdict_register_property(self, PK_BOOL, ptr, name, description);
}


/* String property */
void ptrdict_register_string_property(section_t *self, char *ptr, int maxlen,
				     char *name, char *description)
{
  ptrdict_register_property(self, PK_FORTRAN_STRING, ptr, name,
			   description)->tag = maxlen;
}


/* Point property */
void ptrdict_register_point_property(section_t *self, double *ptr, char *name,
				    char *description)
{
  ptrdict_register_property(self, PK_POINT, ptr, name, description);
}


/* Integer point property */
void ptrdict_register_intpoint_property(section_t *self, int *ptr, char *name,
				       char *description)
{
  ptrdict_register_property(self, PK_INTPOINT, ptr, name, description);
}


/* Enum property */
void ptrdict_register_enum_property(section_t *self, int *ptr, int nchoices,
				   int lenchoice, char *choices, char *name,
				   char *description)
{
  property_t *p;

  p = ptrdict_register_property(self, PK_ENUM, ptr, name, description);

  p->tag = nchoices;
  p->tag2 = lenchoice;
  p->tag4 = choices;
}


/* Variable size list of doubles property */
void ptrdict_register_list_property(section_t *self, double *ptr, int maxlen,
				   int *len, char *name, char *description)
{
  property_t *p;

#ifdef DEBUG
  printf("ptrdict_register_list_property: %p %p\n", self, *self);
#endif

  p = ptrdict_register_property(self, PK_LIST, ptr, name, description);

  p->tag2 = maxlen;
  p->tag5 = len;
}


/* Variable size list of doubles property */
void ptrdict_register_string_list_property(section_t *self, char *ptr,
					  int strlen, int maxlen, int *len,
					  char *name, char *description)
{
  property_t *p;

#ifdef DEBUG
  printf("ptrdict_register_string_list_property: %p %p\n", self, *self);
#endif

  p = ptrdict_register_property(self, PK_FORTRAN_STRING_LIST, ptr, name, description);

  p->tag = strlen;
  p->tag2 = maxlen;
  p->tag5 = len;
}


/* Variable size list of integers property */
void ptrdict_register_integer_list_property(section_t *self, double *ptr,
					   int maxlen, int *len, char *name,
					   char *description)
{
  property_t *p;

#ifdef DEBUG
  printf("ptrdict_register_integer_list_property: %p %p\n", self, *self);
#endif

  p = ptrdict_register_property(self, PK_INT_LIST, ptr, name, description);

  p->tag2 = maxlen;
  p->tag5 = len;
}


/* 1D array property */
void ptrdict_register_array1d_property(section_t *self, double *ptr, int nx,
                                       char *name, char *description)
{
  property_t *p;

  p = ptrdict_register_property(self, PK_ARRAY1D, ptr, name, description);

  p->tag = nx;
}


/* 2D array property */
void ptrdict_register_array2d_property(section_t *self, double *ptr, int nx,
				                               int ny, char *name, char *description)
{
  property_t *p;

  p = ptrdict_register_property(self, PK_ARRAY2D, ptr, name, description);

  p->tag = nx;
  p->tag2 = ny;
}


/* 3D array property */
void ptrdict_register_array3d_property(section_t *self, double *ptr, int nx,
				                               int ny, int nz, char *name,
				                               char *description)
{
  property_t *p;

  p = ptrdict_register_property(self, PK_ARRAY3D, ptr, name, description);

  p->tag = nx;
  p->tag2 = ny;
  p->tag3 = nz;
}


/* Clean up, remove everything from memory */
void ptrdict_cleanup(section_t *root)
{
  section_t *curs, *nexts;
  property_t *curp, *nextp;

  nextp = root->first_property;
  while (nextp) {
    curp = nextp;
    nextp = curp->next;

    free(curp);
  }

  nexts = root->first_child;
  while (nexts) {
    curs = nexts;
    nexts = curs->next;

    ptrdict_cleanup(curs);
  }

  free(root);
}


/* Find the subsection called name. */
section_t *ptrdict_find_section(section_t *self, char *name)
{
  section_t *s;

  s = self->first_child;

  while (s && strcmp(s->name, name) && !(s->nalias > 0 && !strcmp(s->alias, name)))
    s = s->next;

  if (s && s->kind == SK_1TON && s->callback)
    return (section_t *) s->callback(s);
  else
    return s;
}


/* Find the properties called name. */
property_t *ptrdict_find_property(section_t *self, char *name)
{
  property_t *p;

  p = self->first_property;

  while (p && strcmp(p->name, name))
    p = p->next;

  return p;
}


/* Set the property to a value (provided by a string) */
void ptrdict_set_property(property_t *p, char *value)
{
  int i, j;
  double d;
  char *endptr;

  /* Reading "point" data */
  double *point;
  int *intpoint;
  char *fstr;
  char *str;
  char *c1, *c2, *endptr1, *endptr2;

  if (p->provided) {
    printf("[ptrdict_set_property] Error: Property '%s' of section '%s' has "
           "already been set.\n", p->name, p->parent->name);
    exit(1);
  }

  if (!p->ptr) {
    printf("[ptrdict_set_property] Error: Trying to set property '%s' of "
           "section '%s' which has NULL pointer.\n", p->name, p->parent->name);
    exit(1);
  }

#ifdef DEBUG
  printf("ptrdict_set_property: %s = '%s'\n", p->name, value);
#endif

  switch (p->kind) {
  case PK_INT:
    i = strtol(value, &endptr, 10);

    if (endptr != value+strlen(value)) {
      printf("[ptrdict_set_property] Error: Cannot convert '%s' to integer "
             "for property '%s' of section '%s'.\n", value, p->name,
             p->parent->name);
      exit(1);
    }

    *((int*) p->ptr) = i;
    break;
  case PK_DOUBLE:
    d = strtod(value, &endptr);

    if (endptr != value+strlen(value)) {
      printf("[ptrdict_set_property] Error: Cannot convert '%s' to double for property '%s' of section '%s'.\n",
	     value, p->name, p->parent->name);
      exit(1);
    }

    *((double*) p->ptr) = d;
    break;
  case PK_BOOL:
    if (!strcmp(value, "yes") || !strcmp(value, "true") || !strcmp(value, "1")) {
      *((BOOL*) p->ptr) = TRUE;
    } else if (!strcmp(value, "no") || !strcmp(value, "false") || !strcmp(value, "0")) {
      *((BOOL*) p->ptr) = FALSE;
    } else {
      printf("[ptrdict_set_property] Error: Cannot convert '%s' to logical value for property '%s' of section '%s'. "
	     "Valid logical values are 'yes', 'true', '1', 'no', 'false' or '0'.\n",
	     value, p->name, p->parent->name);
      exit(1);
    }
    break;
  case PK_STRING:
    strncpy((char*) p->ptr, value, p->tag-1);
    break;
  case PK_FORTRAN_STRING:
    /* Fortran strings... Fill with blanks. */
    strncpy((char*) p->ptr, value, p->tag);

    if (strlen(value) < p->tag) {
      memset(((char*) p->ptr)+strlen(value), ' ', p->tag-strlen(value));
    } else {
      printf("[ptrdict_set_property] Error: String '%s' in property '%s' of section '%s' exceeds maximum length of %i.\n", value, p->name, p->parent->name, p->tag);
      exit(1);
    };
    break;
  case PK_POINT:
    point = (double*) p->ptr;
    str = strdup(value);
  
    c1 = strchr(str, ',');
    c2 = strrchr(str, ',');

    if (c1 == NULL || c2 == NULL || c1 == c2) {
      printf("[ptrdict_set_property] Error: The point property '%s' of section '%s' needs to have the format 'x, y, z'.\n", p->name, p->parent->name);
      free(str);
      exit(1);
    }

    c1[0] = 0;
    c2[0] = 0;

    while (isblank(*(++c1)));
    while (isblank(*(++c2)));

    point[0] = strtod(str, &endptr);
    point[1] = strtod(c1, &endptr1);
    point[2] = strtod(c2, &endptr2);

    if (endptr != str+strlen(str) || endptr1 != c1+strlen(c1) || endptr2 != c2+strlen(c2)) {
      printf("[ptrdict_set_property] Error: Could not convert property '%s' of section '%s' to a list of numbers.\n", p->name, p->parent->name);
      free(str);
      exit(1);
    }
    free(str);
    break;
  case PK_INTPOINT:
    intpoint = (int*) p->ptr;
    str = strdup(value);
  
    c1 = strchr(str, ',');
    c2 = strrchr(str, ',');

    if (c1 == NULL || c2 == NULL || c1 == c2) {
      printf("[ptrdict_set_property] Error: The point property '%s' of section '%s' needs to have the format 'x, y, z'.\n", p->name, p->parent->name);
      free(str);
      exit(1);
    }

    c1[0] = 0;
    c2[0] = 0;

    while (isblank(*(++c1)));
    while (isblank(*(++c2)));

    intpoint[0] = strtol(str, &endptr, 10);
    intpoint[1] = strtol(c1, &endptr1, 10);
    intpoint[2] = strtol(c2, &endptr2, 10);

    if (endptr != str+strlen(str) || endptr1 != c1+strlen(c1) || endptr2 != c2+strlen(c2)) {
      printf("[ptrdict_set_property] Error: Could not convert property '%s' of section '%s' to a list of numbers.\n", p->name, p->parent->name);
      free(str);
      exit(1);
    }
    break;
  case PK_ENUM:
    j = -1;
    for (i = 0; i < p->tag; i++) {
      if (!strcmp(value, p->tag4 + i*p->tag2))
        j = i;
    };

    if (j < 0) {
      printf("[ptrdict_set_property] Error: Could not find key '%s' in property '%s' of section '%s'.\n",
             value, p->name, p->parent->name);
      exit(1);
    } else {
      *((int*) p->ptr) = j;
    };
    break;
  case PK_LIST:
    point = (double*) p->ptr;
    str = strdup(value);
    c2 = str;

    i = 0;
    c1 = strchr(c2, ',');
    while (c1 && i < p->tag2) {
      c1[0] = 0;
      point[i] = strtod(c2, &endptr);

      while (isblank(*(++c1)));
      c2 = c1;
      c1 = strchr(c2, ',');
      i++;
    }
    if (i >= p->tag2) {
      printf("[ptrdict_set_property] Too many values for property '%s' of "
	     "section '%s'.\n",
             p->name, p->parent->name);
    }
    point[i] = strtod(c2, &endptr);

    *p->tag5 = i+1;

    free(str);
    break;
  case PK_INT_LIST:
    intpoint = (int*) p->ptr;
    str = strdup(value);
    c2 = str;

    i = 0;
    c1 = strchr(c2, ',');
    while (c1 && i < p->tag2) {
      c1[0] = 0;
      intpoint[i] = strtod(c2, &endptr);

      while (isblank(*(++c1)));
      c2 = c1;
      c1 = strchr(c2, ',');
      i++;
    }
    if (i >= p->tag2) {
      printf("[ptrdict_set_property] Too many values for property '%s' of "
	     "section '%s'.\n",
             p->name, p->parent->name);
    }
    intpoint[i] = strtod(c2, &endptr);

    *p->tag5 = i+1;

    free(str);
    break;
  case PK_FORTRAN_STRING_LIST:
    fstr = (char*) p->ptr;

    memset(fstr, ' ', p->tag*p->tag2);

    str = strdup(value);
    c2 = str;

    i = 0;
    c1 = strchr(str, ',');
    while (c1 && i < p->tag2) {
      c1[0] = 0;
      strncpy(fstr + i*p->tag, c2, MIN(p->tag, strlen(c2)));

      while (isblank(*(++c1)));
      c2 = c1;
      c1 = strchr(c2, ',');
      i++;
    }
    if (i >= p->tag2) {
      printf("[ptrdict_set_property] Too many values for property '%s' of section '%s'.\n",
             p->name, p->parent->name);
    }
    strncpy(fstr + i*p->tag, c2, MIN(p->tag, strlen(c2)));

    *p->tag5 = i+1;

    free(str);
    break;
  case PK_ARRAY2D:
    point = (double*) p->ptr;
    str = strdup(value);
    c2 = str;

    i = 0;
    j = p->tag*p->tag2;
    c1 = strchr(c2, ',');
    while (c1 && i < j) {
      c1[0] = 0;
      point[i] = strtod(c2, &endptr);

      while (isblank(*(++c1)));
      c2 = c1;
      c1 = strchr(c2, ',');
      i++;
    }
    if (i >= j) {
      printf("[ptrdict_set_property] Too many values for property '%s' of "
	     "section '%s'.\n",
             p->name, p->parent->name);
    }
    else if (i == 0) {
        /* Fill whole array with this value. */
        point[0] = strtod(c2, &endptr);
        for (i = 1; i < j; i++)  point[i] = point[0];
    }
    else if (i < j-1) {
      printf("[ptrdict_set_property] Too few values for property '%s' of "
             "section '%s'.\n",
             p->name, p->parent->name);
    }
    else {
        point[i] = strtod(c2, &endptr);
    }

    free(str);
    break;
  case PK_ARRAY3D:
    point = (double*) p->ptr;
    str = strdup(value);
    c2 = str;

    i = 0;
    j = p->tag*p->tag2*p->tag3;
    c1 = strchr(c2, ',');
    while (c1 && i < j) {
      c1[0] = 0;
      point[i] = strtod(c2, &endptr);

      while (isblank(*(++c1)));
      c2 = c1;
      c1 = strchr(c2, ',');
      i++;
    }
    if (i >= j) {
      printf("[ptrdict_set_property] Too many values for property '%s' of "
             "section '%s'.\n",
             p->name, p->parent->name);
    }
    else if (i == 0) {
        /* Fill whole array with this value. */
        point[0] = strtod(c2, &endptr);
        for (i = 1; i < j; i++)  point[i] = point[0];
    }
    else if (i < j-1) {
      printf("[ptrdict_set_property] Too few values for property '%s' of "
             "section '%s'.\n",
             p->name, p->parent->name);
    }
    else {
        point[i] = strtod(c2, &endptr);
    }

    free(str);
    break;
  default:
    printf("[ptrdict_set_property] Internal error: Unknown kind %i of property '%s' of section '%s'.\n",
	   p->kind, p->name, p->parent->name);
    exit(1);
  }

  p->provided = TRUE;
}


#define MAX_STR 100

/* Write the property into a string */
void ptrdict_get_property(property_t *p, char *value, int n)
{
  double *point;
  int *intpoint;
  int i, j;
  char str[MAX_STR];

  switch (p->kind) {
  case PK_INT:
    snprintf(value, n, "%i", *((int*) p->ptr));
    break;
  case PK_DOUBLE:
    snprintf(value, n, "%e", *((double*) p->ptr));
    break;
  case PK_BOOL:
    if (*((BOOL*) p->ptr))
      strcpy(value, "yes");
    else
      strcpy(value, "no");
    break;
  case PK_STRING:
    strncpy(value, (char*) p->ptr, MIN(p->tag, n));
    value[MIN(p->tag, n)-1] = 0;
    break;
  case PK_FORTRAN_STRING:
    memcpy(value, (char*) p->ptr, MIN(p->tag, n));
    i = MIN(p->tag, n)-1;
    while (i >= 0 && ((char*) p->ptr)[i] == ' ')  i--;
    value[i+1] = 0;
    break;
  case PK_POINT:
    point = (double*) p->ptr;
    snprintf(value, n, "%f, %f, %f", point[0], point[1], point[2]);
    break;
  case PK_INTPOINT:
    intpoint = (int*) p->ptr;
    snprintf(value, n, "%i, %i, %i", intpoint[0], intpoint[1], intpoint[2]);
    break;
  case PK_ENUM:
    i = *((int*) p->ptr);
    if (i < 0 || i >= p->tag) {
      printf("[ptrdict_get_property] Internal error: Value %i does not map to string for property '%s' of section '%s'.\n",
             i, p->name, p->parent->name);
      exit(1);
    };
    strncpy(value, p->tag4 + i*p->tag2, MIN(p->tag2, n));
    break;
  case PK_LIST:
    point = (double*) p->ptr;
    strcpy(value, "");
    for (i = 0; i < *p->tag5; i++) {
      if (i == 0) 
        snprintf(value, n, "%f", point[i]);
      else
        snprintf(value, n, "%s, %f", value, point[i]);
    }
    break;
  case PK_INT_LIST:
    intpoint = (int*) p->ptr;
    strcpy(value, "");
    for (i = 0; i < *p->tag5; i++) {
      if (i == 0) 
        snprintf(value, n, "%i", intpoint[i]);
      else
        snprintf(value, n, "%s, %i", value, intpoint[i]);
    }
    break;
  case PK_FORTRAN_STRING_LIST:
    strcpy(value, "");
    for (i = 0; i < *p->tag5; i++) {
      memcpy(str, (char*) p->ptr + i*p->tag, MIN(p->tag, MAX_STR));
      j = MIN(p->tag, MAX_STR)-1;
      while (j >= 0 && str[j] == ' ')  j--;
      str[j+1] = 0;
      if (i == 0)
        snprintf(value, n, "%s", str);
      else
        snprintf(value, n, "%s, %s", value, str);
    }
    break;
  default:
    printf("[ptrdict_get_property] Internal error: Unknown kind %i of property '%s' of section '%s'.\n",
	   p->kind, p->name, p->parent->name);
    exit(1);
  }
}


/* Write the names of all subsection to stream f. */
void ptrdict_enum_subsections(section_t *self, FILE *f)
{
  section_t *s;

  if (self->first_child) {
    s = self->first_child;
    while (s) {
      fprintf(f, "- %s\n", s->name);
      fprintf(f, "%s\n", s->description);
      s = s->next;
    }
  } else {
    fprintf(f, "Section '%s' does not have subsections.\n", self->name);
  }
}


/* Write the names of all properties to stream f. */
void ptrdict_enum_properties(section_t *self, FILE *f)
{
  property_t *p;

  if (self->first_property) {
    p = self->first_property;
    while (p) {
      fprintf(f, "- %s\n", p->name);
      fprintf(f, "%s\n", p->description);
      p = p->next;
    }
  } else {
    fprintf(f, "Section '%s' does not have properties.", self->name);
  }
}

/* --- IO ------------------------------------------------------------------------ */

#define iskeywordchar(c)  (isalnum(c) || c == '_')


BOOL skipline(parser_t *p)
{
  while (fgetc(p->f) != '\n' && !feof(p->f));

  p->row++;

  return feof(p->f);
}


char read_char(parser_t *p)
{
  char c;

  if (feof(p->f)) {
    c = '0';
  } else {
    c = fgetc(p->f);
    while (c == '#' || c == '!') {
      skipline(p);
      if (feof(p->f)) {
	c = '0';
      } else {
	c = fgetc(p->f);
      };
    }
  }

  return c;
}


void read_next_token(parser_t *p, char *s, size_t n)
{
  int i = 0;
  char c, last_c;

  while (isblank(c = read_char(p)))
    p->column++;

  last_c = 0;
  while (c != 0 && c != -1 && !isspace(c) && !iskeywordchar(c) && c != '"' && last_c != ';') {
    s[i] = c;
    //    printf("'%c' %i\n", c, c);
    i++;

    if (i >= n) {
      s[i] = 0;
      printf("[read_next_token] Error: Token too long in line %i. Token = '%s'\n", p->row, s);
      exit(1);
    }

    last_c = c;
    c = read_char(p);
    p->column++;
  }
  s[i] = 0;

  if (c != 0 && c != -1) 
    ungetc(c, p->f);
}


void read_next_keyword(parser_t *p, char *s, size_t n, BOOL *is_token)
{
  int i = 0;
  char c;

  *is_token = FALSE;

  while (isblank(c = read_char(p)))
    p->column++;

  if (iskeywordchar(c)) {
  	/* This is really a keyword */
  	while (iskeywordchar(c)) {
      s[i] = c;
      i++;

      if (i >= n) {
        printf("[read_next_keyword] Error: Keyword too long in line %i.\n", p->row);
        exit(1);
      }

      c = read_char(p);
      p->column++;
    }
    s[i] = 0;

    ungetc(c, p->f);
  } else {
  	/* This might be a token */
  	*is_token = TRUE;
  	ungetc(c, p->f);
  	p->column--;
  	
  	read_next_token(p, s, n);
  }
}


void read_next_value(parser_t *p, char *s, size_t n)
{
  int i = 0;
  char c;

  while (isblank(c = read_char(p)))
    p->column++;

  if (c != '"') {
    printf("[read_next_value] Error: \" expected in line %i.\n", p->row);
    exit(1);
  }

  c = fgetc(p->f);
  p->column++;
  while (c != '"' && c != '\n' && c != '\r') {
    s[i] = c;
    i++;

    if (i >= n) {
      s[i] = 0;
      printf("[read_next_value] Error: Token too long in line %i. Token = '%s'\n", p->row, s);
      exit(1);
    }

    c = fgetc(p->f);
    p->column++;
  }
  s[i] = 0;
  p->column++;

  if (c != '"') {
    printf("[read_next_value] Error: End-of-line encountered before '\"' in line %i.\n", p->row);
    exit(1);
  }
}

void finish_line(parser_t *p)
{
  char c;

  while (isblank(c = read_char(p)) && !feof(p->f))
    p->column++;

  if (c != '\n' && !feof(p->f)) {
    printf("[finish_line] Error: End-of-line expected in line %i.\n", p->row);
    exit(1);
  }

  p->column = 1;
  p->row++;
}

int at_end_of_file(parser_t *p)
{
  char c;
  while (isspace(c = read_char(p)) && !feof(p->f)) {
    p->column++;
    if (c == '\n') {
      p->column = 1;
      p->row++;
    }
  }

  if (feof(p->f))
    return 1;
  else {
    ungetc(c, p->f);

    return 0;
  }
}



#define MAX_KEYWORD  130
#define MAX_VALUE  1024*8
#define MAX_TOKEN  1024
#define MAX_SPACE  100
#define INDENT_INC 2


/*
 * These are the long-format parser...
 */

/* Read the current section until 'endsection' is reached. */
void ptrdict_lf_read_section(section_t *self, parser_t *parser)
{
  BOOL section_done = FALSE;
  section_t *s;
  property_t *p;
  char keyword[MAX_KEYWORD+1];
  char value[MAX_VALUE+1];
  char token[MAX_TOKEN+1];
  BOOL is_token;

  self->provided = TRUE;
  if (self->provided_notification)
    *self->provided_notification = TRUE;

  while (!section_done && !at_end_of_file(parser)) {
    read_next_keyword(parser, keyword, MAX_KEYWORD, &is_token);
    
    if (is_token) {
      printf("[ptrdict_lf_read_section] Keyword expected in line %i.\n", parser->row);
      exit(1);
    }

    if (!strcmp(keyword, "endsection")) {
      read_next_value(parser, value, MAX_VALUE);

      if (strcmp(value, self->name)) {
	printf("[ptrdict_lf_read_section] Current open section is '%s', cannot close section '%s' in line %i.\n",
	       self->name, value, parser->row);
	exit(1);
      }

      if (self->kind != SK_SECTION) {
	printf("[ptrdict_lf_read_section] Current open object '%s' is not a section (line %i).\n",
	       self->name, parser->row);
	exit(1);
      }

      finish_line(parser);
      section_done = TRUE;
    } else if (!strcmp(keyword, "endmodule")) {
      read_next_value(parser, value, MAX_VALUE);

      if (strcmp(value, self->name)) {
	printf("[ptrdict_lf_read_section] Current open module is '%s', cannot close module '%s' in line %i.\n",
	       self->name, value, parser->row);
	exit(1);
      }

      if (self->kind != SK_MODULE) {
	printf("[ptrdict_lf_read_section] Current open object '%s' is not a section (line %i).\n",
	       self->name, parser->row);
	exit(1);
      }

      finish_line(parser);
      section_done = TRUE;
    } else if (!strcmp(keyword, "section")) {
      read_next_value(parser, value, MAX_VALUE);
      finish_line(parser);

      s = ptrdict_find_section(self, value);
      if (!s) {
	printf("[ptrdict_lf_read_section] Unknown section '%s' in line %i.\n", value, parser->row);
	
	ptrdict_enum_subsections(self, stdout);

	exit(1);
      } else if (s->kind == SK_MODULE) {
	printf("[ptrdict_lf_read_section] '%s' is a module identifier (line %i).\n", value, parser->row);
	exit(1);
      }

      ptrdict_lf_read_section(s, parser);
    } else if (!strcmp(keyword, "module")) {
      read_next_value(parser, value, MAX_VALUE);
      finish_line(parser);

      s = ptrdict_find_section(self, value);
      if (!s) {
	printf("[ptrdict_lf_read_section] Unknown section '%s' in line %i.\n", value, parser->row);
	
	ptrdict_enum_subsections(self, stdout);

	exit(1);
      } else if (s->kind == SK_SECTION) {
	printf("[ptrdict_lf_read_section] '%s' is a section identifier (line %i).\n", value, parser->row);
	exit(1);
      }

      ptrdict_lf_read_section(s, parser);
    } else {
      p = ptrdict_find_property(self, keyword);

      if (!p) {
	printf("[ptrdict_lf_read_section] Unknown keyword '%s' in line %i.\n"
	       "Possibilities are 'section', 'module' or one of the properties of section '%s', which are:\n",
	       keyword, parser->row, self->name);

	ptrdict_enum_properties(self, stdout);

	exit(1);
      }

      read_next_token(parser, token, MAX_TOKEN);

      if (strcmp(token, "=")) {
	printf("[ptrdict_lf_read_section] '=' expected for assignment of property '%s' in line %i.\n", keyword, parser->row);
	exit(1);
      }

      read_next_value(parser, value, MAX_VALUE);
      finish_line(parser);

      ptrdict_set_property(p, value);
    }
  }

  if (!section_done) {
    printf("[ptrdict_lf_read_section] Error: End-of-file reached, but keyword 'endsection' is missing (line %i).\n", parser->row);
    exit(1);
  }
}


/* Write the current section. */
void ptrdict_write_section(section_t *self, FILE *f, int indent)
{
  int i;
  char space[MAX_SPACE+1], space2[MAX_SPACE+1];
  char value[MAX_VALUE+1];
  section_t *s;
  property_t *p;

#ifdef DEBUG
  printf("ptrdict_write_section: %s, %i\n", self->name, self->kind);
#endif

  if (self->kind == SK_1TON) {

    s = self->first_child;
    while (s) {
      ptrdict_write_section(s, f, indent);

      s = s->next;
    }

  } else if (!(self->kind == SK_MODULE && !self->provided)) {
    if (indent > MAX_SPACE-INDENT_INC) {
      printf("[ptrdict_write_section] Internal error: indent too large.\n");
      exit(1);
    }

    for (i = 0; i < indent+INDENT_INC; ++i) {
      space[i] = ' '; space2[i] = ' ';
    }
    space[indent] = 0;
    space2[indent+INDENT_INC] = 0;

    //fprintf(f, "\n%s# %s\n", space, self->description);

    fprintf(f, "%s%s {\n", space, self->name);

    p = self->first_property;
    while (p) {
      if (p->kind != PK_ARRAY2D && p->kind != PK_ARRAY3D) {
	ptrdict_get_property(p, value, MAX_VALUE);

	fprintf(f, "\n%s  # %s\n", space2, p->description);
	fprintf(f, "%s  %s = \"%s\";\n", space2, p->name, value);
      }

      p = p->next;
    }

    s = self->first_child;
    while (s) {
      ptrdict_write_section(s, f, indent+2*INDENT_INC);

      s = s->next;
    }

    fprintf(f, "%s};\n", space);
  }
}


/* Write ptrdict to a file. */
void ptrdict_write(section_t *root, char *fn)
{
  FILE *f;

  f = fopen(fn, "w");

  ptrdict_write_section(root, f, 0);

  fclose(f);
}



/*
 * These are the short-format parser...
 */

/* Read the current section until 'endsection' is reached. */
void ptrdict_sf_read_section(section_t *self, parser_t *parser)
{
  BOOL section_done = FALSE;
  section_t *s;
  property_t *p;
  char keyword[MAX_KEYWORD+1];
  char value[MAX_VALUE+1];
  char token[MAX_TOKEN+1];
  BOOL is_token;

  self->provided = TRUE;
  if (self->provided_notification)
    *self->provided_notification = TRUE;

  while (!section_done && !at_end_of_file(parser)) {
    read_next_keyword(parser, keyword, MAX_KEYWORD, &is_token);
    if (is_token)
      strcpy(token, keyword);
    else
      read_next_token(parser, token, MAX_TOKEN);
    
    if (!strcmp(token, "};")) {
      section_done = TRUE;
    } else if (!strcmp(token, "}")) {
      read_next_token(parser, token, MAX_TOKEN);
      if (strcmp(token, ";")) {
      	printf("[ptrdict_sf_read_section] ';' expected in line %i.\n", parser->row);
      	exit(1);
      }
      section_done = TRUE;
    } else if (!strcmp(token, "{};")) {
      /* This is an empty module. */
      s = ptrdict_find_section(self, keyword);
      if (!s) {
	printf("[ptrdict_sf_read_section] Unknown section '%s' in line %i.\n", keyword, parser->row);
	
	ptrdict_enum_subsections(self, stdout);

	exit(1);
      }
      
      if (s->kind != SK_MODULE) {
      	printf("[ptrdict_sf_read_section] Module expected, but section encountered in line %i.", parser->row);
      	exit(1);
      }
      
      s->provided = TRUE;
  	  if (s->provided_notification)
    	*s->provided_notification = TRUE;
    } else if (!strcmp(token, "{")) {
      /* We have a section or module */
      s = ptrdict_find_section(self, keyword);
      if (!s) {
	printf("[ptrdict_sf_read_section] Unknown section '%s' in line %i.\n", keyword, parser->row);
	
	ptrdict_enum_subsections(self, stdout);

	exit(1);
      }

      ptrdict_sf_read_section(s, parser);
   } else if (!strcmp(token, "=")) {
      p = ptrdict_find_property(self, keyword);

      if (!p) {
	printf("[ptrdict_sf_read_section] Unknown property '%s' of section '%s' in line %i.\n"
	       "Possibilities are:\n",
	       keyword, self->name, parser->row);
      
	ptrdict_enum_properties(self, stdout);

	exit(1);
      }

      read_next_value(parser, value, MAX_VALUE);

      ptrdict_set_property(p, value);
      
      read_next_token(parser, token, MAX_TOKEN);
      
      if (strcmp(token, ";")) {
      	printf("[ptrdict_sf_read_section] ';' expected in line %i.\n", parser->row);
      	
      	exit(1);
      }
    } else {
      printf("[ptrdict_sf_read_section] Syntax error in line %i. Token = '%s'\n", parser->row, token);
      exit(1);
    }
  }

  if (!section_done) {
    printf("[ptrdict_sf_read_section] Error: End-of-file reached, but file is incomplete.");
    exit(1);
  }
}


/* Read ptrdict from a stream. */
void ptrdict_from_stream(section_t *root, FILE *f)
{
  /* Note: We might want to switch to XML, i.e., using libXML2 eventually.
     This will make this stuff a lot easier, too. */

  parser_t p;
  section_t *s;
  char keyword[MAX_KEYWORD+1];
  char token[MAX_TOKEN+1];
  char value[MAX_VALUE+1];
  BOOL is_token;

  p.column = 1;
  p.row = 1;
  p.f = f;

  /* First line needs to be section with the name of the root section. */
  read_next_keyword(&p, keyword, MAX_KEYWORD, &is_token);
  if (!strcmp(keyword, "section")) {
  	/* Okay, it's the long format. */
  	
  	read_next_value(&p, value, MAX_VALUE);
  	if (strcmp(value, root->name)) {
      printf("[ptrdict_read] Error: Expected '%s' as the name of the first section in line %i.\n", root->name, p.row);
      exit(1);
  	}

    ptrdict_lf_read_section(root, &p);  	
  	
  } else if (!strcmp(keyword, root->name)) {
  	/* Okay, it's the short format. */
  	
  	read_next_token(&p, token, MAX_TOKEN);
  	
  	if (strcmp(token, "{")) {
  	  printf("[ptrdict_read] '{' expected in line %i.\n", p.row);
  	  exit(1);
  	}
  	
	ptrdict_sf_read_section(root, &p);
	  
  } else {
    printf("[ptrdict_read] Error: Keyword 'section' or '%s' expected in line %i.\n", root->name, p.row);
    exit(1);
  }
}


/* Read ptrdict from a file. */
void ptrdict_read(section_t *root, char *fn)
{
  FILE *f = fopen(fn, "r");

  if (!f) {
    printf("[ptrdict_read] Error opening file '%s'.\n%s\n", fn, strerror(errno));
    exit(1);
  }

  ptrdict_from_stream(root, f);

  fclose(f);
}


#ifndef __APPLE__
/* Read ptrdict from a string. */
void ptrdict_from_string(section_t *root, char *s)
{
  FILE *f = fmemopen(s, strlen(s), "r");

  if (f) {
    printf("[ptrdict_read] Something went wrong during fmemopen.\n%s\n",
	   strerror(errno));
    exit(1);
  }

  ptrdict_from_stream(root, f);

  fclose(f);
}
#endif
