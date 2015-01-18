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
#include <stdlib.h>

#include "ptrdict.h"

#include "%(dispatch_header)s"


/*
 * %(disclaimer)s
 */



/*
 * Prototypes
 */

%(prototypes)s



/*
 * Classes
 */

%(classes)s



/*
 * Registration, instantiation
 */

section_t *%(name)s_callback(section_t *hook)
{
  %(name)s_class_t *class;
  section_t *section;
  void *dispatch;

  class = (%(name)s_class_t *) hook->tag;
  dispatch = (%(name)s_class_t *) hook->tag2;

  class->new_instance(dispatch, hook, &section);

  if (strcmp(hook->name, section->name)) {
    printf("[%(name)s_callback] Class and section name differ "
	   "('%%s' != '%%s').\n", hook->name, section->name);
    exit(1);
  }

  return section;
}


void %(name)s_factory_register(section_t *self, void *dispatch)
{
  int i;
  section_t *s;

  for (i = 0; i < N_CLASSES; i++) {
    s = (section_t *)  ptrdict_register_group(self, SK_1TON,
			  		      %(name)s_classes[i].name,
					      "(placeholder)", NULL);
    s->callback = (callback_t) %(name)s_callback;
    s->tag = (void *) &%(name)s_classes[i];
    s->tag2 = (void *) dispatch;
  }
}

