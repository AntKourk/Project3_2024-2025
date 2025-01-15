#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include "inputs.h"
#include "flipEdges.h"
#include "output.h"
#include "centroid.h"
#include "projection.h"
#include "center.h"
#include "circumcenter.h"
#include "inside_convex_polygon_centroid.h"
#include "local_search.h"
#include "simulated_annealing.h"
#include "ant_colony.h"

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;


using namespace std;

bool detectCycle(int node, int parent, const unordered_map<int, vector<int>>& adjacency, unordered_set<int>& visited) {
    visited.insert(node);
    for (int neighbor : adjacency.at(node)) {
        if (neighbor == parent) continue;  // Skip the edge that led here
        if (visited.find(neighbor) != visited.end() || detectCycle(neighbor, node, adjacency, visited)) {
            return true;
        }
    }
    return false;
}

bool hasClosedConstraints(const vector<vector<int>>& shape) {
    unordered_map<int, vector<int>> adjacency;

    // Build adjacency list for constraint edges
    for (const auto& polygon : shape) {
        adjacency[polygon[0]].push_back(polygon[1]);
        adjacency[polygon[1]].push_back(polygon[0]);
    }

    // Detect cycles in the graph using visited set
    unordered_set<int> visited;
    for (const auto& node_edges : adjacency) {
        int node = node_edges.first;
        if (visited.find(node) == visited.end()) {
            if (detectCycle(node, -1, adjacency, visited)) {
                return true;  // Found a cycle -> closed constraint
            }
        }
    }
    return false;  // No cycles -> open constraints
}

int main(int argc, char* argv[]) {
    std::string input_file;
    std::string output_file;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-i" && i + 1 < argc) {
            input_file = argv[++i];
        } else if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else {
            std::cerr << "Usage: " << argv[0] << " -i /path/to/input.json -o /path/to/output.json" << std::endl;
            return 1;
        }
    }

    if (input_file.empty() || output_file.empty()) {
        std::cerr << "Input and output file paths must be provided." << std::endl;
        return 1;
    }

    // Call the inputs function with the input file path
    InputData input = inputs(input_file);

    // Initialize the Constrained Delaunay Triangulation (CDT)
    CDT cdt;

    // Get points
    vector<Point> points = input.points;

    vector<int> points_x = input.points_x;
    vector<int> points_y = input.points_y;


    // Get Region Boundary
    vector<int> region_boundary = input.region_boundary;

    // Get method
    std::string method = input.method;

    // Get parameters
    std::map<std::string, double> parameters = input.parameters;

    double alpha, beta, xi, psi, lambda, kappa, L; 

    // Get if the triangulation is Delaunay or not
    bool isDelaunay = input.delaunay;

    // Insert points into the triangulation
    for (const Point& p : points) {
        cdt.insert(p);
    }

    // Insert the region boundary as a constrained polygon
    std::vector<Point> polygon;
    for (int idx : region_boundary) {
        if (idx < points.size()) {
            polygon.push_back(points[idx]);
        } else {
            cerr << "Invalid index in region_boundary: " << idx << endl;
        }
    }

    // Append the first point again to close the polygon
    if (!region_boundary.empty()) {
        int first_idx = region_boundary[0];
        if (first_idx < points.size()) {
            polygon.push_back(points[first_idx]);
        } else {
            cerr << "Invalid first index in region_boundary: " << first_idx << endl;
        }
    }

    // Check if the polygon is valid and insert the constraint
    if (polygon.size() > 2) {
        cdt.insert_constraint(polygon.begin(), polygon.end());
    } else {
        cerr << "Not enough points to form a boundary." << endl;
    }

    // Define and add the constrained edges (from additional_constraints)
    const std::vector<std::vector<int>>& constraints = input.additional_constraints;

    // Insert constrained edges based on the provided indices
    for (const auto& constraint : constraints) {  
        if (constraint.size() == 2) {
            int idx1 = constraint[0];
            int idx2 = constraint[1];
            if (idx1 < points.size() && idx2 < points.size()) {
                cdt.insert_constraint(points[idx1], points[idx2]);
            } else {
                cerr << "Invalid constraint index: " << idx1 << ", " << idx2 << endl;
            }
        }
    }

    std::vector<std::vector<int>> shape;

    // Add region boundary edges as pairs
    for (size_t i = 0; i < region_boundary.size(); ++i) {
        int start = region_boundary[i];
        int end = region_boundary[(i + 1) % region_boundary.size()];  
        shape.push_back({start, end});
    }

    // Both region boundry and constraints
    shape.insert(shape.end(), constraints.begin(), constraints.end());

    char category;

    Polygon_2 is_convex_polygon;
    for (const auto& pt : polygon) {
        is_convex_polygon.push_back(pt);
    }

    if (is_convex_polygon.is_convex()) {
        if(constraints.size() == 0){
            category = 'A';
        } else {
            category = !hasClosedConstraints(shape) ? 'B' : 'C';
        }
    } else {
        category = 'E';
        if(constraints.size() == 0){
            bool flag = false;
            for(int i=0; i<region_boundary.size(); i++) {
                int current = region_boundary[i];
                int next = region_boundary[(i + 1) % region_boundary.size()];
                if (points_x[current] == points_x[next] || points_y[current] == points_y[next]) {
                    flag = true;
                    break;
                }
            }
            if(flag == true){
                category = 'D';
            }
        }
    }

    cout << "Category: " << category << endl; 

    if(!isDelaunay){
        // Prompt user to choose the Steiner point insertion method
        cout << "Please choose a method for Steiner points from the following options:\n";
        cout << "1: Center of longest edge\n";
        cout << "2: Projection\n";
        cout << "3: Circumcenter\n";
        cout << "4: Centroid of internal convex polygon\n";
        cout << "5: Centroid\n";
        cout << "6: Flip\n";
        cout << "Enter the number corresponding to your choice: ";

        int choice;
        cin >> choice;

        // Execute the chosen method based on user input
        switch (choice) {
            case 1:
                cdt = center_steiner_points(points, cdt);
                break;
            case 2:
                cdt = projection(points, cdt);
                break;
            case 3:
                cdt = circumcenter_steiner_points(points, cdt);
                break;
            case 4:
                cdt = inside_convex_polygon_centroid_steiner_points(points, cdt);
                break;
            case 5:
                cdt = centroid_steiner_points(points, cdt);
                break;
            case 6:
                cdt = flip_edges(points, cdt);
                break;
            default:
                cerr << "Invalid choice. Please enter a number between 1 and 5.\n";
                return 1;
        }
        if (method == "local") {
            L = parameters["L"];
            local_search(points, cdt, L, input_file, output_file);
        } else if (method == "sa") {
            alpha = parameters["alpha"];
            beta = parameters["beta"];
            L = parameters["L"];
            simulated_annealing(points, cdt, alpha, beta, L, input_file, output_file);
        } else if (method == "ant") {
            alpha = parameters["alpha"];
            beta = parameters["beta"];
            xi = parameters["xi"];
            psi = parameters["psi"];
            lambda = parameters["lambda"];
            kappa = parameters["kappa"];
            L = parameters["L"];
            ant_colony(points, cdt, L, kappa, alpha, beta, lambda, xi, psi, input_file, output_file);
        } else {
            cerr << "Invalid method specified." << endl;
            return 1;
        }
    
    } else {
        if (method == "local") {
            L = parameters["L"];
            local_search(points, cdt, L, input_file, output_file);
        } else if (method == "sa") {
            alpha = parameters["alpha"];
            beta = parameters["beta"];
            L = parameters["L"];
            simulated_annealing(points, cdt, alpha, beta, L, input_file, output_file);
        } else if (method == "ant") {
            alpha = parameters["alpha"];
            beta = parameters["beta"];
            xi = parameters["xi"];
            psi = parameters["psi"];
            lambda = parameters["lambda"];
            kappa = parameters["kappa"];
            L = parameters["L"];
            ant_colony(points, cdt, L, kappa, alpha, beta, lambda, xi, psi, input_file, output_file);
        } else {
            cerr << "Invalid method specified." << endl;
            return 1;
        }
    }

    return 0;
}
