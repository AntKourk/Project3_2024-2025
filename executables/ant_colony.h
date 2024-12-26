#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>

// Define CGAL kernel and types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

int ant_colony(std::vector<Point> points, DT& dt, int L, int kappa, double alpha, double beta, double lamda, double xi, double psi, const std::string& input_file, const std::string& output_file);
K::FT height(const Point& P, const Point& A, const Point& B);
K::FT radius_to_height_ratio(const FaceHandle& face, DT dt);
double calculateDelta(DT& dt, double alpha, double beta, int steiner_points_count);
std::vector<double> steiner_point_probability(const std::vector<double>& t, const std::vector<double>& h, double x, double y);
double calculateEnergyAnt(DT& dt, double alpha, double beta, int steiner_points_count);
int steiner_method(const std::vector<double>& probabilities);