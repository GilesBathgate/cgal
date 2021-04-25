// Copyright (c) 1997-2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Courtesy of LEDA
#ifndef CGAL_CHAINED_MAP_H
#define CGAL_CHAINED_MAP_H

#include <CGAL/memory.h>
#include <iostream>

namespace CGAL {

namespace internal {

template <typename T, typename Allocator = CGAL_ALLOCATOR(T) >
class chained_map
{
   struct element_type
   {
     std::size_t key;
     T value;
     element_type*  next;
   };

public:
   static constexpr std::size_t min_size = 32;
   typedef element_type* item;

   T& inf(item it) const { return it->value; }
   T& xdef() { return stop.value; }
   const T& cxdef() const { return stop.value; }

   chained_map(std::size_t n = min_size)
   {
       init_table(n);
   }

   void clear()
   {
       for(item p = table_begin + 1; p < table_free; p++)
         if (p->key != nullptr_key || p >= table_begin + table_size)
           p->value = T();

       for (item item = table_begin ; item != table_end ; ++item)
           destroy(item);
       alloc.deallocate(table_begin, table_end - table_begin);

       init_table(min_size);
   }

   ~chained_map()
   {
     for (item item = table_begin ; item != table_end ; ++item)
       destroy(item);
     alloc.deallocate(table_begin, table_end - table_begin);
   }

   T& access(item it, std::size_t key)
   {
       stop.key = key;
       item q = it->next;
       while (q->key != key) q = q->next;
       if (q != &stop)
       { old_index = key;
           return q->value;
       }

       // index x not present, insert it

       if (table_free == table_end)   // table full: rehash
       { rehash();
           it = hash(key);
       }

       if (it->key == nullptr_key)
       { it->key = key;
           it->value = stop.value;  // initializes p->i to xdef
           return it->value;
       }

       q = table_free++;
       q->key = key;
       q->value = stop.value;    // initializes q->i to xdef
       q->next = it->next;
       it->next = q;
       return q->value;
   }

   T& access(std::size_t key)
   { item p = hash(key);

       if ( p->key == key ) {
           old_index = key;
           return p->value;
       }
       else {
           if ( p->key == nullptr_key ) {
               p->key = key;
               p->value = stop.value;  // initializes p->i to xdef
               old_index = key;
               return p->value;
           } else
               return access(p,key);
       }
   }

   item lookup(std::size_t key) const
   { item p = hash(key);
       ((std::size_t &)stop.key) = key;  // cast away const
       while (p->key != key)
       { p = p->next; }
       return (p == &stop) ? 0 : p;
   }

   void statistics() const
   { std::cout << "table_size: " << table_size <<"\n";
       std::size_t n = 0;
       for (item p = table_begin + 1; p < table_begin + table_size; p++)
           if (p ->key != nullptr_key) n++;
       std::size_t used_in_overflow = table_free - (table_begin + table_size );
       n += used_in_overflow;
       std::cout << "number of entries: " << n << "\n";
       std::cout << "fraction of entries in first position: " <<
                    ((double) (n - used_in_overflow))/n <<"\n";
       std::cout << "fraction of empty lists: " <<
                    ((double) (n - used_in_overflow))/table_size<<"\n";
   }

private:

   typedef std::allocator_traits<Allocator> allocator_traits;
   typedef typename allocator_traits::template rebind_alloc<element_type> allocator_type;
   typedef std::allocator_traits<allocator_type> allocator_type_traits;

   item  hash(std::size_t key)  const
   { return table_begin + (key & table_size_1);  }

   void init_table(std::size_t n)
   {
       std::size_t t = min_size;
       while (t < n) t <<= 1;

       table_size = t;
       table_size_1 = t-1;
       std::size_t s = t + t/2;
       table_begin = alloc.allocate(s);
       for (std::size_t i = 0 ; i < s ; ++i){
           allocator_type_traits::construct(alloc,table_begin + i);
       }

       table_free = table_begin + t;
       table_end = table_begin + s;

       for (item p = table_begin; p < table_free; p++)
       { p->next = &stop;
           p->key = nullptr_key;
       }
       table_begin->key = non_nullptr_key;
   }

   void rehash()
   {
       item old_table_begin = table_begin;
       item old_table_end = table_end;
       std::size_t old_table_size = table_size;
       std::size_t old_table_size_1 = table_size_1;
       item old_table_free = table_free;

       item old_table_mid = table_begin + table_size;

       init_table(2*table_size);

       item p;

       for(p = old_table_begin + 1; p < old_table_mid; p++)
       { std::size_t x = p->key;
           if ( x != nullptr_key ) // list p is non-empty
           { item q = hash(x);
               q->key = x;
               q->value = p->value;
           }
       }

       while (p < old_table_end)
       { std::size_t x = p->key;
           insert(x,p->value);
           p++;
       }

       // delete old table
       item save_table = table_begin;
       item save_table_end = table_end;
       item save_free = table_free;
       std::size_t save_table_size = table_size;
       std::size_t save_table_size_1 = table_size_1;

       table_begin = old_table_begin;
       table_end = old_table_end;
       table_size = old_table_size;
       table_size_1 = old_table_size_1;
       table_free = old_table_free;

       old_table_begin = 0;

       T v = access(old_index);

       for (item item = table_begin ; item != table_end ; ++item)
           destroy(item);
       alloc.deallocate(table_begin, table_end - table_begin);

       table_begin = save_table;
       table_end = save_table_end;
       table_size = save_table_size;
       table_size_1 = save_table_size_1;
       table_free = save_free;
       access(old_index) = v;
   }

   inline void insert(std::size_t key, T value)
   { item q = hash(key);
       if ( q->key == nullptr_key ) {
           q->key = key;
           q->value = value;
       } else {
           table_free->key = key;
           table_free->value = value;
           table_free->next = q->next;
           q->next = table_free++;
       }
   }

   void destroy(item item)
   {
     allocator_type_traits::destroy(alloc,item);
   }

   static constexpr std::size_t nullptr_key = 0;
   static constexpr std::size_t non_nullptr_key = 1;

   element_type stop;
   item table_begin;
   item table_end;
   item table_free;
   std::size_t table_size;
   std::size_t table_size_1;
   std::size_t old_index;
   allocator_type alloc;
};

} // namespace internal
} //namespace CGAL

#endif // CGAL_CHAINED_MAP_H
