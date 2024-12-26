#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include "output.h"
#include <string>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;
typedef K::FT FT;
typedef CGAL::Polygon_2<K> Polygon;

//Function to check whether a point is inside convex hull
bool is_within_convex_hull(const Point& p, const std::vector<Point>& convex_hull) {
    Polygon hull_polygon(convex_hull.begin(), convex_hull.end());
    return CGAL::bounded_side_2(hull_polygon.vertices_begin(), hull_polygon.vertices_end(), p, K()) 
           != CGAL::ON_UNBOUNDED_SIDE;
}

// Function to calculate the angle between two points and a common vertex
template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3) {
    double a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3)));
    double c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

    // Cosine rule to calculate the angle
    double cos_angle = (b * b + c * c - a * a) / (2 * b * c);
    return std::acos(cos_angle) * 180.0 / M_PI; // Return angle in degrees
}

// Function to check if a triangle is obtuse
template <typename FaceHandle>
int obtuse_vertex_index(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    
    if (angle1 > 90.0 + 0.01) return 0;
    if (angle2 > 90.0 + 0.01) return 1;
    if (angle3 > 90.0 + 0.01) return 2;
    return -1;
}

// Function to calculate the circumcenter of a triangle
Point circumcenter(const Point& p1, const Point& p2, const Point& p3) {
    // Using the CGAL function to calculate the circumcenter
    return CGAL::circumcenter(p1, p2, p3);
}

template <typename DT>
std::pair<std::vector<Point>, std::vector<Point>> add_steiner_in_centroid(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points,  const std::vector<Point>& convex_hull) {
    bool added_steiner = false;

    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            Point p1 = face->vertex(0)->point();
            Point p2 = face->vertex(1)->point();
            Point p3 = face->vertex(2)->point();
            Point circumcenter_point = circumcenter(p1, p2, p3);

            if (is_within_convex_hull(circumcenter_point, convex_hull)) {
                steiner_points.push_back(circumcenter_point);
                points.push_back(circumcenter_point);
                added_steiner = true;
            }
        }
    }

    // Insert Steiner points into the triangulation
    for (const Point& p : steiner_points) {
        dt.insert(p);
    }

    return {steiner_points, points};
}

DT circumcenter_steiner_points(std::vector<Point> points, DT dt) {

    std::vector<Point> convex_hull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull));
    std::vector<std::pair<size_t, size_t>> edges;

    bool obtuse_exists = true;
    int iterations = 0;
    int obtuse_counter=0, previous_obtuse_counter=0, counter=0;
    
    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    while (obtuse_exists && iterations <= 5 && counter < 3) {
        all_points = add_steiner_in_centroid(dt, steiner_points, points, convex_hull);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        
        obtuse_exists = false;

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            int obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
                obtuse_counter++;
            }
        }
        iterations++;
        if(obtuse_counter > previous_obtuse_counter){
            counter++;
            previous_obtuse_counter = obtuse_counter;
        }
    }

    return dt;
}
