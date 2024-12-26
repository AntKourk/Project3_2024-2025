#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <ctime>
#include <cmath>    // For pow
#include <random>   // For random number generation
#include "ant_colony.h"
#include "projection.h"
#include "circumcenter.h"
#include "center.h"
#include "inside_convex_polygon_centroid.h" 
#include "utils.h"
#include "output.h"
#include <string>

// Define CGAL kernel and types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

// Function to compute the height of point P above line segment AB
K::FT height(const Point& P, const Point& A, const Point& B) {
    // Calculate the determinant for the area of the parallelogram
    K::FT det = (B.x() - A.x()) * (A.y() - P.y()) - (A.x() - P.x()) * (B.y() - A.y());
    
    // Length of the line segment AB
    K::FT length_squared = (B.x() - A.x()) * (B.x() - A.x()) + (B.y() - A.y()) * (B.y() - A.y());

    if (length_squared == K::FT(0)) {
        throw std::runtime_error("Segment length is zero; cannot compute height.");
    }

    double length = CGAL::to_double(length_squared);
    double abs_det = CGAL::to_double(CGAL::abs(det));

    // Return the height as a K::FT
    return K::FT(abs_det / std::sqrt(length));
}

// Function to calculate the radius-to-height ratio
K::FT radius_to_height_ratio(const FaceHandle& face, DT dt) {
    Point A = face->vertex(0)->point();
    Point B = face->vertex(1)->point();
    Point C = face->vertex(2)->point();

    int obtuse_adjacent_triangles_counter = 0;

    if (obtuse_vertex_index(face) != -1) {
            obtuse_adjacent_triangles_counter++;
            // Explore neighboring faces
            for (int i = 0; i < 3; ++i) {
                FaceHandle neighbor_face = face->neighbor(i);
                if (!dt.is_infinite(neighbor_face) && obtuse_vertex_index(neighbor_face) != -1) {
                    obtuse_adjacent_triangles_counter++;
                }
            }
        }

    K::FT area = CGAL::abs(CGAL::area(A, B, C));
    if (area == K::FT(0)) {
        throw std::runtime_error("Degenerate triangle; cannot compute circumradius.");
    }

    K::FT circumradius_squared = CGAL::squared_distance(A, B) * CGAL::squared_distance(A, C) * CGAL::squared_distance(B, C)
                                  / (4 * area);

    double circumradius_double = std::sqrt(CGAL::to_double(circumradius_squared));

    K::FT circumradius = K::FT(circumradius_double);

    K::FT h = height(C, A, B);
    if (h == K::FT(0)) {
        throw std::runtime_error("Height is zero; invalid triangle.");
    }

    K::FT r = circumradius / h;

    //heuristics based on steiner point techniques
    if (obtuse_adjacent_triangles_counter >=2) {
        return K::FT(1);
    } else if (r >= 1 && r <= 2) {
        return r / (2.0 + r);
    } else if (r < 1) {
        return (3.0 - 2.0*r) / 3.0 > K::FT(0) ? (3.0 - 2.0*r) / 3.0 : K::FT(0);
    } else {
        return (r - 1.0) / r > K::FT(0) ? (r - 1.0) / r : K::FT(0);
    }
}

double calculateDelta(DT& dt, double alpha, double beta, int steiner_points_count) {
    int obtuse_count = 0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++obtuse_count;
        }
    }

    return 1.0/(1 + alpha * obtuse_count + beta * steiner_points_count);
}

// Function to calculate probabilities
std::vector<double> steiner_point_probability(const std::vector<double>& t, const std::vector<double>& h, double x, double y, std::vector<double> probabilities) {

    std::vector<double> prob(t.size());

    // Compute weighted values for each option
    double denominator = 0.0;
    for (size_t i = 0; i < t.size(); i++) {
        //probabilities[i] = std::pow(t[i], x) * std::pow(h[i], y);
        denominator += probabilities[i];
    }

    // Normalize probabilities
    for (size_t i = 0; i < t.size(); i++) {
        prob[i] = probabilities[i]/denominator;
    }

    return prob;
}


// Function to select an option based on probabilities
int steiner_method(const std::vector<double>& probabilities) {
    // Compute the sum of probabilities
    double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    
    // Normalize the probabilities if the sum is not 0
    std::vector<double> normalized_probabilities(probabilities.size());
    for (size_t i = 0; i < probabilities.size(); ++i) {
        normalized_probabilities[i] = probabilities[i] / sum;
    }

    // Compute cumulative probabilities
    std::vector<double> cumulative(normalized_probabilities.size());
    std::partial_sum(normalized_probabilities.begin(), normalized_probabilities.end(), cumulative.begin());

    // Generate a random number in the range [0, 1)
    std::random_device rd;  // Seed generator
    std::mt19937 gen(rd()); // Random number engine
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double random_value = dis(gen);

    // Select the option based on cumulative probabilities
    for (size_t i = 0; i < cumulative.size(); i++) {
        if (random_value < cumulative[i]) {
            return i; // Return the selected option index
        }
    }

    // Fallback (should not happen if probabilities are valid)
    return -1;
}

double calculateEnergyAnt(DT& dt, double alpha, double beta, int steiner_points_count) {
    int obtuse_count = 0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++obtuse_count;
        }
    }

    if(obtuse_count == 0){
       return 0.0; 
    }
    else {
        return alpha * obtuse_count + beta * steiner_points_count;
    }
    
}

// Ant colony optimization function
int ant_colony(std::vector<Point> points, DT& dt, int L, int kappa, double alpha, double beta, double lamda, double xi, double psi, const std::string& input_file, const std::string& output_file ) {
    bool obtuse_exists = true;
    int iterations = 0;
    int obtuse_count = 0, obtuse_previous_count = 0;

    std::vector<double> t(4, 1.0);      // Initialize with size 4, default value 1.0
    std::vector<double> deltaT(4, 0.0); // Initialize with size 4, default value 0.0
    std::vector<double> h(4, 1.0);
    std::vector<double> probabilities(4, 1.0);
    std::vector<double> total_probabilities(4, 1.0);

    DT temp_dt;

    int method_used = 0;

    double sum = 0.0, previous_energy, new_energy, deltaE;

    std::vector<Point> steiner_points;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt); // Draw initial triangulation

    previous_energy = 100;
    
    obtuse_exists = false;

    for (int cycle = 1; cycle < L; cycle++) {
        obtuse_previous_count = 0;
        obtuse_count = 0;
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_previous_count++;
            }
        }

        for (int ant = 1; ant < kappa; ant++) {
            for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
                int obtuse_vertex = obtuse_vertex_index(face);
                if (obtuse_vertex != -1) {
                    Point p_obtuse = face->vertex(obtuse_vertex)->point();
                    Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
                    Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();

                    total_probabilities = steiner_point_probability(t, h, xi, psi,probabilities);

                    method_used = steiner_method(probabilities);

                    Point new_point;
                    switch (method_used) {
                        case 0:
                            new_point = project_point_onto_line(p_obtuse, p1, p2);
                            break;
                        case 1:
                            new_point = circumcenter(p_obtuse, p1, p2);
                            break;
                        case 2:
                            new_point = longest_edge_center(p1, p2);
                            break;
                        case 3:
                            auto polygon_points = find_convex_polygon(dt, face);
                            if (!polygon_points.empty()) {  // Only proceed if convex polygon found
                                new_point = compute_centroid(polygon_points);
                            }
                            break;
                    }

                    //Calculate the energy's reduction
                    temp_dt.insert(new_point);
                    h[method_used] = CGAL::to_double(radius_to_height_ratio(face, dt));
                    probabilities[method_used] = std::pow(t[method_used], xi) * std::pow(h[method_used], psi);
                    new_energy = calculateEnergyAnt(temp_dt, alpha, beta, steiner_points.size()+1);
                    deltaE = new_energy - previous_energy;

                    if(deltaE < 0) {
                        previous_energy = new_energy;
                        points.push_back(new_point);
                        steiner_points.push_back(new_point);
                        dt.insert(new_point);
                    }

                    temp_dt = dt;
                    break;
                }
            }   
        }

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_count++;
            }
        }

        if(obtuse_count < obtuse_previous_count){
            deltaT[method_used] = calculateDelta(dt, alpha, beta, steiner_points.size());
        } else {
            deltaT[method_used] = 0.0;
        }
        t[method_used] = (1 - lamda) * t[method_used] + deltaT[method_used];
    }

    edges = print_edges(dt, points);
    output(edges, steiner_points, input_file, output_file, obtuse_previous_count);
    CGAL::draw(dt); // Draw final triangulation

    return 0;
}
