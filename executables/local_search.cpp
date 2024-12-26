#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
#include "output.h"
#include "utils.h"
#include "local_search.h"
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

// Find best Steiner
std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points) {
    bool added_steiner = false;

        // Store obtuse triangles
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            int obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                // Get the vertices of the obtuse triangle
                Point p_obtuse = face->vertex(obtuse_vertex)->point();
                Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
                Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();

                // Test multiple methods for Steiner point placement
                Point projection_point = project_point_onto_line(p_obtuse, p1, p2);
                Point circumcenter_point = circumcenter(p_obtuse, p1, p2);
                Point centroid_point = calculate_centroid(p_obtuse, p1, p2);
                Point center = longest_edge_center(p1, p2);
                Point inside_convex_polygon_centroid;
                auto polygon_points = find_convex_polygon(dt, face);
                    if (!polygon_points.empty()) {  // Only proceed if convex polygon found
                        inside_convex_polygon_centroid = compute_centroid(polygon_points);
                    }

                // Evaluate obtuse triangle reduction for each method
                DT temp_dt = dt;
                temp_dt.insert(projection_point);
                int count_projection = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(circumcenter_point);
                int count_circumcenter = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(centroid_point);
                int count_centroid = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(inside_convex_polygon_centroid);
                int count_inside_convex_polygon_centroid = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(center);
                int count_center = count_obtuse_triangles(temp_dt);

                // Select the best point
                Point best_point;
                int min_count = std::min({count_projection, count_circumcenter, count_centroid, count_center, count_inside_convex_polygon_centroid});
                if (min_count == count_projection) {
                    best_point = projection_point;
                } else if (min_count == count_circumcenter) {
                    best_point = circumcenter_point;
                } else if (min_count == count_center) {
                    best_point = center;
                } else if (min_count == count_inside_convex_polygon_centroid) {
                    best_point = inside_convex_polygon_centroid;
                } else {
                    best_point = centroid_point;
                }

                // Add the best point to the DT
                points.push_back(best_point);
                steiner_points.push_back(best_point);
                added_steiner = true;
            }
        }

    for (const Point& p : steiner_points) {
        dt.insert(p);
    }

    return {steiner_points, points};
}

int local_search(std::vector<Point> points, DT dt, int max_iterations, const std::string& input_file, const std::string& output_file) {
    bool obtuse_exists = true;
    int obtuse_count = 0, obtuse_previous_count = 0;
    int iterations = 0;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    int counter = 0;
    while (obtuse_exists && iterations <= max_iterations && counter < 3) {
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_previous_count++;
            }
        }
        all_points = add_best_steiner(dt, steiner_points, points);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        obtuse_exists = false;
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
                obtuse_count++;
            }
        }
        iterations++;
        if (obtuse_count > obtuse_previous_count) {
            counter++;
        }
    }

    edges = print_edges(dt, all_points.first);
    output(edges, steiner_points, input_file, output_file, obtuse_count);
    CGAL::draw(dt);

    return 0;
}