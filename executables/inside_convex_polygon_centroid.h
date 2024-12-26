#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;
typedef DT::Face_handle FaceHandle;

Point compute_centroid(const std::vector<Point>& points);

std::vector<Point> find_convex_polygon(DT& dt, FaceHandle start_face);

DT inside_convex_polygon_centroid_steiner_points(std::vector<Point> points, DT dt);