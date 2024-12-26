#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <cmath>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;


double calculateEnergy(DT& dt, double alpha, double beta, int steiner_points_count);
int generate_random_number();
bool accept_new_configuration(double deltaE, double T);
// Function to perform Simulated Annealing
int simulated_annealing(std::vector<Point> points, DT dt, double alpha, double beta, int L, const std::string& input_file, const std::string& output_file);

