#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/bind.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "output.h"
#include "inputs.h"
#include <string>
#include <sstream>
#include <boost/algorithm/string/replace.hpp>
#include <fstream>
#include <map>
#include <set>
#include <gmpxx.h>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;

void write_json_no_escaping(const boost::property_tree::ptree& pt, const std::string& filename) {
    std::ostringstream oss;
    write_json(oss, pt, true); 

    std::string json_str = oss.str();
    
    // Replace the escaped `\/` with `/`
    boost::replace_all(json_str, "\\/", "/");

    // Write the modified JSON string to the file
    std::ofstream outfile(filename);
    if (outfile) {
        outfile << json_str;
        outfile.close();
        std::cout << "Output written to " << filename << std::endl;
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
}

// Function to convert K::FT to a string representation
std::string rational_to_string(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    if (exact_coord.get_den() == 1) {
        return exact_coord.get_num().get_str();
    } else {
        std::ostringstream oss;
        oss << exact_coord.get_num().get_str() << "/" << exact_coord.get_den().get_str();
        return oss.str();
    }
}

// Function to convert a map into a boost::property_tree::ptree
boost::property_tree::ptree map_to_ptree(const std::map<std::string, double>& map) {
    boost::property_tree::ptree pt;
    for (const auto& [key, value] : map) {
        pt.put(key, value);
    }
    return pt;
}

void output(const std::vector<std::pair<size_t, size_t>>& edges, std::vector<Point> steiner_points_given, const std::string& input_file, const std::string& output_file, int obtuse_count) {
    // Creation of property tree
    boost::property_tree::ptree pt;

    // Read from JSON file
    try {
        read_json(input_file, pt);
    } catch (const boost::property_tree::json_parser_error &e) {
        std::cerr << "Error reading JSON: " << e.what() << std::endl;
    }

    // Get data from property tree
    std::string instance_uid = pt.get<std::string>("instance_uid");
    std::string method = pt.get<std::string>("method");

    // Specify known integer parameters
    std::set<std::string> int_parameters = {"kappa", "L"};

    // Populate the parameters map with either integer or double based on the known types
    std::map<std::string, double> parameters;
    for (const auto& item : pt.get_child("parameters")) {
        if (int_parameters.find(item.first) != int_parameters.end()) {
            parameters[item.first] = static_cast<double>(item.second.get_value<int>());
        } else {
            parameters[item.first] = item.second.get_value<double>();
        }
    }

    std::vector<Point> steiner_points = steiner_points_given;

    // Create the output property tree
    boost::property_tree::ptree output_pt;

    // Populate the JSON structure
    output_pt.put("content_type", "CG_SHOP_2025_Solution");
    output_pt.put("instance_uid", instance_uid);

    // Steiner points x
    boost::property_tree::ptree steiner_points_x_node;
    for (const auto& x : steiner_points) {
        boost::property_tree::ptree child;
        child.put("", rational_to_string(x.x()));
        steiner_points_x_node.push_back(std::make_pair("", child));
    }
    output_pt.add_child("steiner_points_x", steiner_points_x_node);

    // Steiner points y
    boost::property_tree::ptree steiner_points_y_node;
    for (const auto& y : steiner_points) {
        boost::property_tree::ptree child;
        child.put("", rational_to_string(y.y()));
        steiner_points_y_node.push_back(std::make_pair("", child));
    }
    output_pt.add_child("steiner_points_y", steiner_points_y_node);

    // Edges
    boost::property_tree::ptree edges_node;
    for (const auto& edge : edges) {
        boost::property_tree::ptree edge_array;

        // Add the indices to the array
        edge_array.push_back(boost::property_tree::ptree::value_type("", boost::property_tree::ptree(std::to_string(edge.first))));
        edge_array.push_back(boost::property_tree::ptree::value_type("", boost::property_tree::ptree(std::to_string(edge.second))));

        // Add the array to edges_node
        edges_node.push_back(std::make_pair("", edge_array));
    }

    // Add edges_node to the main output property tree
    output_pt.add_child("edges", edges_node);

    output_pt.put("obtuse_count", obtuse_count);

    output_pt.put("method", method);

    // Add the parameters as a child
    output_pt.add_child("parameters", map_to_ptree(parameters));

    // Write the output JSON to a file
    try {
        write_json_no_escaping(output_pt, output_file);
    } catch (const boost::property_tree::json_parser_error &e) {
        std::cerr << "Error writing JSON: " << e.what() << std::endl;
    }
}
