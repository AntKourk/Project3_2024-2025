#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
// #include "output.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;

// Function to calculate the angle between two points and a common vertex
template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3) {
    double a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3)));
    double c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

    // Cosine rule to calculate the angle
    double cos_angle = (b * b + c * c - a * a) / (2 * b * c);
    return std::acos(cos_angle) * 180.0 / M_PI;
}

// Function to check if a triangle is obtuse (has an angle > 90 degrees)
template <typename FaceHandle>
bool is_obtuse_triangle(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    
    // Check if any angle is greater than 90 degrees
    return (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0);
}

// Function to flip the diagonal if there are obtuse triangles
template <typename DT>
void flip_if_obtuse(DT& dt) {
    for (auto edge = dt.finite_edges_begin(); edge != dt.finite_edges_end(); ++edge) {
        FaceHandle face1 = edge->first;
        FaceHandle face2 = face1->neighbor(edge->second);
        
        // Check if any of the triangles is obtuse
        if (is_obtuse_triangle(face1) || is_obtuse_triangle(face2)) {
            // Try flipping the diagonal
            if (dt.is_flipable(edge->first, edge->second)) {
                dt.flip(edge->first, edge->second);
            }
        }
    }
}

DT flip_edges(std::vector<Point> points, DT dt) {
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    // Flip obtuse edges if possible
    flip_if_obtuse(dt);

    return dt;
}