#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;  // This defines the Point type as CGAL's Point_2 type

struct InputData {
    std::vector<Point> points;
    std::vector<int> region_boundary;
    std::vector<std::vector<int>> additional_constraints;
    std::string method;
    std::map<std::string, double> parameters;
    bool delaunay;
};

InputData inputs(const std::string& input_file); 
