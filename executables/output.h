#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <string>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;

void output(DT dt, const std::vector<std::pair<size_t, size_t>>& edges, std::vector<Point> steiner_points_given, const std::string& input_file, const std::string& output_file, int obtuse_count);
