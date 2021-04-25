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
   struct chained_map_elem
   {
     std::size_t k;
     T i;
     chained_map_elem*  succ;
   };

public:
   static constexpr std::size_t min_size = 32;
   typedef chained_map_elem*  chained_map_item;
   typedef chained_map_item item;

   T& inf(chained_map_item it) const { return it->i; }
   T& xdef() { return STOP.i; }
   const T& cxdef() const { return STOP.i; }

   chained_map(std::size_t n = 1) :
       nullptrKEY(0), NONnullptrKEY(1), old_table(0)
   {
       if (n < min_size)
           init_table(min_size);
       else {
           std::size_t ts = 1;
           while (ts < n) ts <<= 1;
           init_table(ts);
       }
   }

   void clear()
   {
       for(chained_map_item p = table + 1; p < free; p++)
         if (p->k != nullptrKEY || p >= table + table_size)
           p->i = T();

       for (chained_map_item item = table ; item != table_end ; ++item)
           destroy(item);
       alloc.deallocate(table, table_end - table);

       init_table(min_size);
   }

   ~chained_map()
   {
     if (old_table)
     {
       for (chained_map_item item = old_table ; item != old_table_end ; ++item)
         destroy(item);
       alloc.deallocate(old_table, old_table_end - old_table);
     }
     for (chained_map_item item = table ; item != table_end ; ++item)
       destroy(item);
     alloc.deallocate(table, table_end - table);
   }

   T& access(chained_map_item p, std::size_t x)
   {
       STOP.k = x;
       chained_map_item q = p->succ;
       while (q->k != x) q = q->succ;
       if (q != &STOP)
       { old_index = x;
           return q->i;
       }

       // index x not present, insert it

       if (free == table_end)   // table full: rehash
       { rehash();
           p = HASH(x);
       }

       if (p->k == nullptrKEY)
       { p->k = x;
           p->i = STOP.i;  // initializes p->i to xdef
           return p->i;
       }

       q = free++;
       q->k = x;
       q->i = STOP.i;    // initializes q->i to xdef
       q->succ = p->succ;
       p->succ = q;
       return q->i;
   }

   T& access(std::size_t x)
   { chained_map_item p = HASH(x);

       if (old_table) del_old_table();
       if ( p->k == x ) {
           old_index = x;
           return p->i;
       }
       else {
           if ( p->k == nullptrKEY ) {
               p->k = x;
               p->i = STOP.i;  // initializes p->i to xdef
               old_index = x;
               return p->i;
           } else
               return access(p,x);
       }
   }

   chained_map_item lookup(std::size_t x) const
   { chained_map_item p = HASH(x);
       ((std::size_t &)STOP.k) = x;  // cast away const
       while (p->k != x)
       { p = p->succ; }
       return (p == &STOP) ? 0 : p;
   }

   void statistics() const
   { std::cout << "table_size: " << table_size <<"\n";
       std::size_t n = 0;
       for (chained_map_item p = table + 1; p < table + table_size; p++)
           if (p ->k != nullptrKEY) n++;
       std::size_t used_in_overflow = free - (table + table_size );
       n += used_in_overflow;
       std::cout << "number of entries: " << n << "\n";
       std::cout << "fraction of entries in first position: " <<
                    ((double) (n - used_in_overflow))/n <<"\n";
       std::cout << "fraction of empty lists: " <<
                    ((double) (n - used_in_overflow))/table_size<<"\n";
   }

private:

   chained_map_elem*  HASH(std::size_t x)  const
   { return table + (x & table_size_1);  }

   void init_table(std::size_t t)
   {
       table_size = t;
       table_size_1 = t-1;
       table = alloc.allocate(t + t/2);
       for (std::size_t i = 0 ; i < t + t/2 ; ++i){
           std::allocator_traits<allocator_type>::construct(alloc,table + i);
       }

       free = table + t;
       table_end = table + t + t/2;

       for (chained_map_item p = table; p < free; p++)
       { p->succ = &STOP;
           p->k = nullptrKEY;
       }
       table->k = NONnullptrKEY;
   }

   void rehash()
   {
       old_table = table;
       old_table_end = table_end;
       old_table_size = table_size;
       old_table_size_1 = table_size_1;
       old_free = free;

       chained_map_item old_table_mid = table + table_size;

       init_table(2*table_size);

       chained_map_item p;

       for(p = old_table + 1; p < old_table_mid; p++)
       { std::size_t x = p->k;
           if ( x != nullptrKEY ) // list p is non-empty
           { chained_map_item q = HASH(x);
               q->k = x;
               q->i = p->i;
           }
       }

       while (p < old_table_end)
       { std::size_t x = p->k;
           insert(x,p->i);
           p++;
       }
   }

   void del_old_table()
   {
       chained_map_item save_table = table;
       chained_map_item save_table_end = table_end;
       chained_map_item save_free = free;
       std::size_t save_table_size = table_size;
       std::size_t save_table_size_1 = table_size_1;

       table = old_table;
       table_end = old_table_end;
       table_size = old_table_size;
       table_size_1 = old_table_size_1;
       free = old_free;

       old_table = 0;

       T p = access(old_index);

       for (chained_map_item item = table ; item != table_end ; ++item)
           destroy(item);
       alloc.deallocate(table, table_end - table);

       table = save_table;
       table_end = save_table_end;
       table_size = save_table_size;
       table_size_1 = save_table_size_1;
       free = save_free;
       access(old_index) = p;
   }

   inline void insert(std::size_t x, T y)
   { chained_map_item q = HASH(x);
       if ( q->k == nullptrKEY ) {
           q->k = x;
           q->i = y;
       } else {
           free->k = x;
           free->i = y;
           free->succ = q->succ;
           q->succ = free++;
       }
   }

   void destroy(chained_map_elem* item)
   {
     typedef std::allocator_traits<allocator_type> Allocator_type_traits;
     Allocator_type_traits::destroy(alloc,item);
   }

   const std::size_t nullptrKEY;
   const std::size_t NONnullptrKEY;

   chained_map_elem STOP;

   chained_map_elem* table;
   chained_map_elem* table_end;
   chained_map_elem* free;
   std::size_t table_size;
   std::size_t table_size_1;

   chained_map_elem* old_table;
   chained_map_elem* old_table_end;
   chained_map_elem* old_free;
   std::size_t old_table_size;
   std::size_t old_table_size_1;

   std::size_t old_index;
   typedef std::allocator_traits<Allocator> Allocator_traits;
   typedef typename Allocator_traits::template rebind_alloc<chained_map_elem> allocator_type;

   allocator_type alloc;
};

} // namespace internal
} //namespace CGAL

#endif // CGAL_CHAINED_MAP_H
