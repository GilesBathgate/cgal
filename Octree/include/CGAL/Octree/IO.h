
#ifndef OCTREE_IO_H
#define OCTREE_IO_H

#include <CGAL/Octree.h>

#include <iostream>

using std::ostream;

template<class Kernel, class PointRange>
void writeToStream(ostream &os, const CGAL::Octree_node<Kernel, PointRange> &node, int depth = 0) {

  // Set indentation
  for (int i = 0; i < depth; ++i) {
    os << ". ";
  }

  // Print out this node
  os << "Node at coordinates (" << node.location().x() << ", " << node.location().y() << ", " << node.location().z()
     << ") "
     << "containing " << node.num_points() << " points \n";

  // Print out this node's children
  if (!node.is_leaf()) {
    for (int i = 0; i < 8; ++i) {
      writeToStream(os, node.children()[i], depth + 1);
    }
  }
}

template<class Kernel, class PointRange>
ostream &operator<<(ostream &os, const CGAL::Octree_node<Kernel, PointRange> &node) {

  //os << "[ Center: (" << node.location().x() << ", " << node.location().y() << ", " << node.location().z() << ") ]";
  writeToStream(os, node);

  return os;
}

#endif //OCTREE_IO_H
