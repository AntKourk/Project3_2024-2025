#ifndef UTILS_H
#define UTILS_H

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <vector>
#include <map>
#include <vector>
#include <set>
#include <algorithm> // for std::min and std::max

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3);

template <typename FaceHandle>
int obtuse_vertex_index(const FaceHandle& face);

int count_obtuse_triangles(const DT& dt);

template <typename DT>
std::vector<std::pair<size_t, size_t>> print_edges(const DT& dt, std::vector<Point> points) {
    std::vector<std::pair<size_t, size_t>> edges;
    std::set<std::pair<size_t, size_t>> unique_edges;

    // Create a map from Point to its index in the points vector
    std::map<typename DT::Point, size_t> point_index_map;
    for (size_t i = 0; i < points.size(); ++i) {
        point_index_map[points[i]] = i;
    }

    // Print edges and their indices
    for (auto edge = dt.finite_edges_begin(); edge != dt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3)->point();
        auto v2 = edge->first->vertex((edge->second + 2) % 3)->point();
        size_t idx1 = point_index_map[v1];
        size_t idx2 = point_index_map[v2];

        // Normalize the edge
        std::pair<size_t, size_t> normalized_edge = {std::min(idx1, idx2), std::max(idx1, idx2)};

        // Add to the set if it's not already present
        if (unique_edges.insert(normalized_edge).second) {
            edges.push_back(normalized_edge);
        }
    }

    return edges;
}

#endif // UTILS_H
