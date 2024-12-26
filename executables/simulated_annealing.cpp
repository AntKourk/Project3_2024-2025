#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include "output.h"
#include "utils.h"
#include "simulated_annealing.h"
#include "projection.h"
#include "circumcenter.h" 
#include "centroid.h"
#include "center.h"
#include "inside_convex_polygon_centroid.h" 
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;

double calculateEnergy(DT& dt, double alpha, double beta, int steiner_points_count) {
    int obtuse_count = 0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++obtuse_count;
        }
    }

    return alpha * obtuse_count + beta * steiner_points_count;
}

int generate_random_number() {
    static std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_int_distribution<int> dist(1, 5); // Range [1, 5]
    return dist(rng); // Generate and return the random number
}

bool accept_new_configuration(double deltaE, double T) {
    // Initialize random number generator
    static std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution [0, 1]

    if (deltaE < 0) {
        // Accept new configuration unconditionally
        return true;
    } else {
        // Calculate probability e^(-∆E/T)
        double probability = std::exp(-deltaE / T);
        // Generate a random number in [0, 1]
        double random = dist(rng);
        // Accept with probability
        return random < probability;
    }
}

int simulated_annealing(std::vector<Point> points, DT dt, double alpha, double beta, int L, const std::string& input_file, const std::string& output_file) {
    double T = 1.0, previous_energy, new_energy, deltaE;
    int random_number;
    int obtuse_count;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    previous_energy = calculateEnergy(dt, alpha, beta, steiner_points.size());

    while (T >= 0) {

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {

                Point p_obtuse = face->vertex(obtuse_vertex)->point();
                Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
                Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();

                // Generate a random Steiner point
                random_number = generate_random_number();
                Point new_point;
                switch (random_number) {
                    case 1:
                        new_point = project_point_onto_line(p_obtuse, p1, p2);
                        break;
                    case 2:
                        new_point = circumcenter(p_obtuse, p1, p2);
                        break;
                    case 3:
                        new_point = calculate_centroid(p_obtuse, p1, p2);
                        break;
                    case 4:
                        new_point = longest_edge_center(p1, p2);
                        break;
                    case 5:
                        auto polygon_points = find_convex_polygon(dt, face);
                        if (!polygon_points.empty()) {  // Only proceed if convex polygon found
                            new_point = compute_centroid(polygon_points);
                        }
                        break;
                }

                //Calculate the energy's reduction
                new_energy = calculateEnergy(dt, alpha, beta, steiner_points.size());
                deltaE = new_energy - previous_energy;

                if(accept_new_configuration(deltaE, T)) {
                    previous_energy = new_energy;
                    points.push_back(new_point);
                    steiner_points.push_back(new_point);
                    dt.insert(new_point);
                }
            }
        }
        T = T - 1.0/L;
    }

    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        auto obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            obtuse_count++;
        }
    }

    edges = print_edges(dt, points);
    output(edges, steiner_points, input_file, output_file, obtuse_count);
    CGAL::draw(dt);

    return 0;
}