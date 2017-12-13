#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <math.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3  Vector_3;
typedef Kernel::Triangle_3  Triangle_3;
typedef Kernel::Aff_transformation_3  Aff_transformation_3;
typedef Kernel::FT Scalar;
typedef Polyhedron::Facet_const_iterator FacetIterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator HalffacetCirculator;

void exportAsciiSTL(Polyhedron poly)
{

        std::cout << "solid stl\n";

        for(FacetIterator fi=poly.facets_begin(); fi!=poly.facets_end(); ++fi) {
                HalffacetCirculator hc=fi->facet_begin();
                CGAL_assertion(circulator_size(hc)==3);
                Triangle_3 t((hc++)->vertex()->point(),
                             (hc++)->vertex()->point(),
                             (hc++)->vertex()->point());

                Vector_3 n=t.supporting_plane().orthogonal_vector();
                Scalar ls=n.squared_length();
                Scalar l=Scalar(sqrt(CGAL::to_double(ls)));
                Vector_3 un=n/l;

                float x1=CGAL::to_double(t[0].x());
                float y1=CGAL::to_double(t[0].y());
                float z1=CGAL::to_double(t[0].z());
                float x2=CGAL::to_double(t[2].x());
                float y2=CGAL::to_double(t[2].y());
                float z2=CGAL::to_double(t[2].z());
                float x3=CGAL::to_double(t[1].x());
                float y3=CGAL::to_double(t[1].y());
                float z3=CGAL::to_double(t[1].z());
                float nx=CGAL::to_double(un.x());
                float ny=CGAL::to_double(un.y());
                float nz=CGAL::to_double(un.z());

                std::cout << "  facet normal " << nx << " " << ny << " " << nz << "\n";
                std::cout << "    outer loop\n";
                std::cout << "      vertex " << x1 << " " << y1 << " " << z1 << "\n";
                std::cout << "      vertex " << x2 << " " << y2 << " " << z2 << "\n";
                std::cout << "      vertex " << x3 << " " << y3 << " " << z3 << "\n";
                std::cout << "    endloop\n";
                std::cout << "  endfacet\n";
        }

        std::cout << "endsolid stl\n";
}

int main() {
  Polyhedron P;
  std::cin >> P;
  if(P.is_closed()) {
    Nef_polyhedron N1(P);

    Polyhedron P2;
    if(N1.is_simple()) {
      for(int i=0; i<1000; ++i) {
        CGAL::convert_nef_polyhedron_to_polygon_mesh(N1,P2);
      }
      exportAsciiSTL(P2);
    }
    else
      std::cerr << "N1 is not a 2-manifold." << std::endl;
  }
}
