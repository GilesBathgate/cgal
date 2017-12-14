#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <fstream>
#include <iostream>
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

int main(int argc, char** argv)
{
	int N=10000;

	if(argc>1)
		N=atoi(argv[1]);

	// create two nested cubes
	Polyhedron poly_in;
	std::ifstream input("cube.off");
	input >> poly_in;

	const Nef_polyhedron nef(poly_in);

	CGAL::Timer timer;
	timer.start();
	for(int i=0; i<N; ++i) {
		Polyhedron poly_out;
		CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, poly_out);
		if(poly_out.size_of_vertices()!=poly_in.size_of_vertices())
			std::cerr <<"ERROR1\n";
	}
	timer.stop();
	std::cout << "convert_nef_polyhedron_to_polygon_mesh " << timer.time() << "\n";


	timer.reset();
	timer.start();
	for(int i=0; i<N; ++i) {
		Polyhedron poly_out;
		CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, poly_out);
		CGAL::Polygon_mesh_processing::triangulate_faces(poly_out);
		if(poly_out.size_of_vertices()!=poly_in.size_of_vertices())
			std::cerr <<"ERROR2\n";
	}
	timer.stop();
	std::cout << "convert_nef_polyhedron_to_polygon_mesh + triangulate_faces " << timer.time() << "\n";


	timer.reset();
	timer.start();
	for(int i=0; i<N; ++i) {
		Polyhedron poly_out;
		CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, poly_out,true);
		if(poly_out.size_of_vertices()!=poly_in.size_of_vertices())
			std::cerr <<"ERROR3\n";
	}
	timer.stop();
	std::cout << "convert_nef_polyhedron_to_polygon_mesh(true) " << timer.time() << "\n";


	timer.reset();
	timer.start();
	for(int i=0; i<N; ++i) {
		Polyhedron poly_out;
		nef.convert_to_polyhedron(poly_out);
		if(poly_out.size_of_vertices()!=poly_in.size_of_vertices())
			std::cerr <<"ERROR4\n";
	}
	timer.stop();
	std::cout << "convert_to_polyhedron " << timer.time() << "\n";


	return 0;

}
